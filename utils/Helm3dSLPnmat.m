function [An,Anx,Any,Anz] = Helm3dSLPnmat(t,s,k)
% Helmholtz SLPn. see https://github.com/bobbielf2/ZetaTrap3D_Unified/blob/main/potentials/Helm3dSLPmat.m
%

if nargin == 0, test_Helm3dSLPn; return; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
An = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddottn./(r.*rr), s.w*(-1/4/pi));  % including src quadr wei
if nargout>1
  tmp = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r)./(r.*rr), s.w*(-1/4/pi));
  Anx = d1.*tmp; 
  Any = d2.*tmp; 
  Anz = d3.*tmp; 
%   AlapS = Lap3dSLPmat(t,s);
%   [~,AlapSnx] = Lap3dSLPnmat(t,s);
%   test = 1i*imag(Anx) - AlapS.*(k*sin(k*r)./r.*d1) + AlapSnx.*cos(k*r);
end
end

function test_Helm3dSLPn
% not an actual test of helmholtz SLPn
% see if singularity is similar to 2d ...

% order
p = 8; % try <= 14, before condition number kicks in
ref = 1; % 1, 2, 3, 4...
omega = 10;

% a patch on a torus
a = 1; b = 0.5; 
phi = outfun(a,b);
chart = @(t,p) [b*sin(t); (a+b*cos(t)).*sin(p); (a+b*cos(t)).*cos(p)]; 
orig = [0;0]; tpansiz1 = 1/2/ref; tpansiz2 = 1/4/ref;
[x0, w0, D] = gauss(p); [x1, x2] = meshgrid(x0); 
r = chart(orig(1)+tpansiz1*x1(:)',orig(2)+tpansiz2*x2(:)'); 
rc = mean(r,2); r = r - rc;
s = get_high_order_quad(r,p);
sf = get_high_order_quad(r,2*p); % upsampled 

% define close targets (num 27419) & density function
tmp = sqrt(sum(s.w))*(rand(3,100000)-1/2);
t.x = mean(s.x,2) + [1.25/3*tmp(1:2,:);1.25*tmp(3,:)]; % random targets
t.x = t.x + rc; t.x = t.x(:,phi(t.x(1,:),t.x(2,:),t.x(3,:))) - rc;
% rng(2)
% t.nx = rand(size(t.x)); t.nx = t.nx./sqrt(sum(t.nx.^2));
t.nx = zeros(size(t.x)); 
t.nx(1,:) = 1;
fmu = @(x,y,z) sin(3*x)+cos(2*y)+sin(z+1/2)+x.^2+y.^3+z.^4+exp(x.*y.*z);

%
[A,Ax,Ay,Az] = Helm3dSLPnmat(t,s,omega);
Af = Helm3dSLPnmat(t,sf,omega);

%
d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');

%
d1f = bsxfun(@minus,t.x(1,:)',sf.x(1,:)); % 3 coords of displacement matrix (M*N)
d2f = bsxfun(@minus,t.x(2,:)',sf.x(2,:));
d3f = bsxfun(@minus,t.x(3,:)',sf.x(3,:));
rrf = d1f.^2+d2f.^2+d3f.^2;   % dist^2 mat
rf  = sqrt(rrf);
ddottnf = bsxfun(@times,d1f,t.nx(1,:)')+bsxfun(@times,d2f,t.nx(2,:)')+bsxfun(@times,d3f,t.nx(3,:)');

% final check before kernel split quadrature
%%% check 1: smooth part
Ai = -1/(4*pi)*(sin(omega*r) - omega*r.*cos(omega*r)).*ddottn./(r.*rr).*s.w;
Aif = -1/(4*pi)*(sin(omega*rf) - omega*rf.*cos(omega*rf)).*ddottnf./(rf.*rrf).*sf.w;
ui = Ai*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uif = Aif*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
erri = abs(ui-uif)/max(abs(uif(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(erri),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(erri) + " ")

%%% check 2: real part, 1./r coefficient
Ar1 = 1/(4*pi)*omega*sin(omega*r)./r.*ddottn.*s.w;
Arf1 = 1/(4*pi)*omega*sin(omega*rf)./rf.*ddottnf.*sf.w;
u1 = Ar1*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uf1 = Arf1*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
err1 = abs(u1-uf1)/max(abs(uf1(:)));
figure(2),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err1),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err1) + " ")

%%% check 3: real part, ddottn./(r.*rr) coefficient
Ar2 = 1/(4*pi)*cos(omega*r).*s.w;
Arf2 = 1/(4*pi)*cos(omega*rf).*sf.w;
u2 = Ar2*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uf2 = Arf2*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
err2 = abs(u2-uf2)/max(abs(uf2(:)));
figure(3),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err2),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err2) + " ")

keyboard

end

function s = get_high_order_quad(r,orderf)

% assume tensor product grid
order = sqrt(numel(r(1,:)));

% interpolation stuff
[x1 x2] = meshgrid(gauss(order)); 
uvs = [x1(:)';x2(:)'];
[xf,wf,D] = gauss(orderf);
[xf1 xf2] = meshgrid(xf);
uvsf = [xf1(:)';xf2(:)'];
[~,V,R] = interpmat_2d(uvsf,uvs);
wwf = wf(:)*wf(:)'; wwf = wwf(:)'; 

% input of grad SLP
sxf = ((V'\R')'*r')';
tp_ind = reshape(1:orderf^2,[orderf orderf]);
xts = nan(1,orderf^2); yts = xts; zts = xts; % partial t (small circle)
for j=1:orderf
  xk = sxf(1,tp_ind(:,j)); yk = sxf(2,tp_ind(:,j)); zk = sxf(3,tp_ind(:,j));
  xts(tp_ind(:,j)) = D*xk(:); yts(tp_ind(:,j)) = D*yk(:); zts(tp_ind(:,j)) = D*zk(:);
end
rts = [xts(:),yts(:),zts(:)]';
xps = nan(1,orderf^2); yps = xps; zps = xps; % partial p (big circle)
for j=1:orderf
  r = sxf(1,tp_ind(j,:)); y = sxf(2,tp_ind(j,:)); z = sxf(3,tp_ind(j,:));
  xps(tp_ind(j,:)) = D*r(:); yps(tp_ind(j,:)) = D*y(:); zps(tp_ind(j,:)) = D*z(:);
end
rps = [xps(:),yps(:),zps(:)]';
nx =  cross(rps, rts);       % outward normal
sp = sqrt(sum(nx.^2,1)); % speeds
nx = nx./sp;
w = sp.*wwf;

s.x = sxf;
s.nx = nx;
s.w = w;

end

function phi = outfun(a,b)
% define indicator function for torus
my_eps = 0.01; % if want to make this smaller, then need to refine reference soln
phi = @(x,y,z) ((z-a*cos(atan2(y,z))).^2 + (y-a*sin(atan2(y,z))).^2 + x.^2) > b^2+my_eps; 
end