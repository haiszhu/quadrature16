function A = Helm3dDLPmat(t,s,k)
% Helmholtz DLP. see https://github.com/bobbielf2/ZetaTrap3D_Unified/blob/main/potentials/Helm3dDLPmat.m
%

if nargin == 0, test_Helm3dDLP; return; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));
A = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddotsn./(r.*rr), s.w*(1/4/pi));  % including src quadr wei

end

function test_Helm3dDLP
% not an actual test of helmholtz DLP
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
fmu = @(x,y,z) sin(3*x)+cos(2*y)+sin(z+1/2)+x.^2+y.^3+z.^4+exp(x.*y.*z);

%
A = Helm3dDLPmat(t,s,omega);
Af = Helm3dDLPmat(t,sf,omega);

%
d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));

%
d1f = bsxfun(@minus,t.x(1,:)',sf.x(1,:)); % 3 coords of displacement matrix (M*N)
d2f = bsxfun(@minus,t.x(2,:)',sf.x(2,:));
d3f = bsxfun(@minus,t.x(3,:)',sf.x(3,:));
rrf = d1f.^2+d2f.^2+d3f.^2;   % dist^2 mat
rf  = sqrt(rrf);
ddotsnf = bsxfun(@times,d1f,sf.nx(1,:))+bsxfun(@times,d2f,sf.nx(2,:))+bsxfun(@times,d3f,sf.nx(3,:));


% final check before kernel split quadrature
% check 1: smooth part
Ai = 1/(4*pi)*(sin(omega*r) - omega*r.*cos(omega*r)).*ddotsn./(r.*rr).*s.w;
Aif = 1/(4*pi)*(sin(omega*rf) - omega*rf.*cos(omega*rf)).*ddotsnf./(rf.*rrf).*sf.w;
ui = Ai*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uif = Aif*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
erri = abs(ui-uif)/max(abs(uif(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(erri),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(erri) + " ")
% check 2: real part, 1./r coefficient
Ar1 = 1/(4*pi)*omega*sin(omega*r)./r.*ddotsn.*s.w;
Arf1 = 1/(4*pi)*omega*sin(omega*rf)./rf.*ddotsnf.*sf.w;
u1 = Ar1*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uf1 = Arf1*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
err1 = abs(u1-uf1)/max(abs(uf1(:)));
figure(2),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err1),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err1) + " ")
% check 3: real part, ddotsn./(r.*rr) coefficient
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
% check 4: together
Atest = 1i*Ai + 1./r.*Ar1 + ddotsn./(r.*rr).*Ar2;
u = A*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
utest = Atest*fmu(s.x(1,:),s.x(2,:),s.x(3,:))'; 
err = abs(u-utest)/max(abs(utest(:)));
figure(4),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err) + " ")
Aftest = 1i*Aif + 1./rf.*Arf1 + ddotsnf./(rf.*rrf).*Arf2;
uf = Af*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
uftest = Aftest*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))'; 
errf = abs(uf-uftest)/max(abs(uftest(:)));
figure(5),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(errf),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(errf) + " ")
% check 5: 
% A5 = imag(A);
% Af5 = imag(Af);
A5 = A - 1/(4*pi)*(1+(omega*r).^2/2+(omega*r).^4*(1/24-1/6)+(omega*r).^6*(-1/720+1/120)+(omega*r).^8*(1/40320-1/5040)).*ddotsn./(r.*rr).*s.w;
Af5 = Af - 1/(4*pi)*(1+(omega*rf).^2/2+(omega*rf).^4*(1/24-1/6)+(omega*rf).^6*(-1/720+1/120)+(omega*rf).^8*(1/40320-1/5040)).*ddotsnf./(rf.*rrf).*sf.w;
% A5 = 1/(4*pi)*(cos(omega*r)+omega*sin(omega*r).*r-1-(omega*r).^2/2).*ddotsn./(r.*rr).*s.w;
% Af5 = 1/(4*pi)*(cos(omega*rf)+omega*sin(omega*rf).*rf-1-(omega*rf).^2/2).*ddotsnf./(rf.*rrf).*sf.w;
% A5 = 1/(4*pi)*(sin(omega*r).*(omega*r)-(omega*r).^2).*ddotsn./(r.*rr).*s.w;
% Af5 = 1/(4*pi)*(sin(omega*rf).*(omega*rf)-(omega*rf).^2).*ddotsnf./(rf.*rrf).*sf.w;
% A5 = 1/(4*pi)*((omega*r).^4).*ddotsn./(r.*rr).*s.w;
% Af5 = 1/(4*pi)*((omega*rf).^4).*ddotsnf./(rf.*rrf).*sf.w;
% A5 = A - 1/(4*pi)*(cos(omega*r)+omega*sin(omega*r).*r).*ddotsn./(r.*rr).*s.w;
% Af5 = Af - 1/(4*pi)*(cos(omega*rf)+omega*sin(omega*rf).*rf).*ddotsnf./(rf.*rrf).*sf.w;
u5 = A5*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uf5 = Af5*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
err5 = abs(u5-uf5)/max(abs(uf5(:)));
figure(6),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err5),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err5) + " ")

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