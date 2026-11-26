function A  = Helm3dSLPmat(t,s,k)
% Helmholtz SLP.

if nargin == 0, test_Helm3dSLP; return; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
r = sqrt(d1.^2+d2.^2+d3.^2);   % dist mat
A = bsxfun(@times, exp(1i*k*r)./r, s.w*(1/4/pi));  % including src quadr wei

end

function test_Helm3dSLP
% not an actual test of helmholtz SLP
% see if singularity is similar to 2d ...

% order
p = 8; % try <= 14, before condition number kicks in
ref = 1; % 1, 2, 3, 4...

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
A = Helm3dSLPmat(t,s,10);
Af = Helm3dSLPmat(t,sf,10);

% 
u = imag(A)*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uf = imag(Af)*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
% error
err = abs(u-uf)/max(abs(uf(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(err) + " ")

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