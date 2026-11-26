function [An,Anx,Any,Anz] = Helm3dDLPnmat(t,s,k)
% Helmholtz DLP. see https://github.com/bobbielf2/ZetaTrap3D_Unified/blob/main/potentials/Helm3dDLPmat.m
%

if nargin == 0, test_Helm3dDLPn; return; end

d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));
% A = bsxfun(@times, (1-1i*k*r).*exp(1i*k*r).*ddotsn./(r.*rr), s.w*(1/4/pi));  % including src quadr wei
% targ deriv wanted
ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
tndotsn = bsxfun(@times,s.nx(1,:),t.nx(1,:)')+bsxfun(@times,s.nx(2,:),t.nx(2,:)')+bsxfun(@times,s.nx(3,:),t.nx(3,:)');
An = bsxfun(@times, exp(1i*k*r)./(r.*rr).*((k^2-3*(1-1i*k*r)./rr).*ddotsn.*ddottn+(1-1i*k*r).*tndotsn), s.w*(1/4/pi));  % dipole deriv, incl src quad wei
if nargout>1
  Anx = bsxfun(@times, exp(1i*k*r)./(r.*rr).*((k^2-3*(1-1i*k*r)./rr).*ddotsn.*d1+(1-1i*k*r).*s.nx(1,:)), s.w*(1/4/pi)); 
  Any = bsxfun(@times, exp(1i*k*r)./(r.*rr).*((k^2-3*(1-1i*k*r)./rr).*ddotsn.*d2+(1-1i*k*r).*s.nx(2,:)), s.w*(1/4/pi));
  Anz = bsxfun(@times, exp(1i*k*r)./(r.*rr).*((k^2-3*(1-1i*k*r)./rr).*ddotsn.*d3+(1-1i*k*r).*s.nx(3,:)), s.w*(1/4/pi));
end
end

function test_Helm3dDLPn
% not an actual test of helmholtz DLPn
% check singularity instead

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
t.nx = rand(size(t.x)); t.nx = t.nx./sqrt(sum(t.nx.^2));
fmu = @(x,y,z) sin(3*x)+cos(2*y)+sin(z+1/2)+x.^2+y.^3+z.^4+exp(x.*y.*z);

%%% 1st check numerical derivative... for a few targets far away
t0.x = s.x(:,:) + s.nx(:,:);
t0.nx = zeros(3,1); t0.nx(2,:) = 2; % change this to 1,2,3 to check all components
delta = 0.000000001;
t0_delta.x = t0.x + delta*t0.nx; % delta for numerical derivative
A0 = Helm3dDLPmat(t0,s,omega);
A0_delta = Helm3dDLPmat(t0_delta,s,omega);
test = (A0_delta-A0)/delta; % is the sign ok?
test2 = Helm3dDLPnmat(t0,s,omega);
diff = abs(test-test2);
max(diff(:)) % numerical differentiation error
% done

%%% break it down
d1 = bsxfun(@minus,t.x(1,:)',s.x(1,:)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,t.x(2,:)',s.x(2,:));
d3 = bsxfun(@minus,t.x(3,:)',s.x(3,:));
rr = d1.^2+d2.^2+d3.^2;   % dist^2 mat
r  = sqrt(rr);
ddotsn = bsxfun(@times,d1,s.nx(1,:))+bsxfun(@times,d2,s.nx(2,:))+bsxfun(@times,d3,s.nx(3,:));
ddottn = bsxfun(@times,d1,t.nx(1,:)')+bsxfun(@times,d2,t.nx(2,:)')+bsxfun(@times,d3,t.nx(3,:)');
tndotsn = bsxfun(@times,s.nx(1,:),t.nx(1,:)')+bsxfun(@times,s.nx(2,:),t.nx(2,:)')+bsxfun(@times,s.nx(3,:),t.nx(3,:)');
[An,Anx,Any,Anz] = Helm3dDLPnmat(t,s,omega);
my1 = (1-1i*omega*r); my2 = exp(1i*omega*r); my3 = ddotsn./(r.*rr); % break DLP down into 3 pieces, A = Helm3dDLPmat(t,s,omega); A2 = my1.*my2.*my3.*(s.w*(1/4/pi));
my1_dx = -1i*omega*d1./r; my2_dx = 1i*omega*my2.*d1./r; my3_dx = 1./(r.*rr).*(-3*d1./rr.*ddotsn + s.nx(1,:));
Anx_test = my1_dx.*my2.*my3.*(s.w*(1/4/pi)) + my1.*my2_dx.*my3.*(s.w*(1/4/pi)) + my1.*my2.*my3_dx.*(s.w*(1/4/pi));
my1_dy = -1i*omega*d2./r; my2_dy = 1i*omega*my2.*d2./r; my3_dy = 1./(r.*rr).*(-3*d2./rr.*ddotsn + s.nx(2,:));
Any_test = my1_dy.*my2.*my3.*(s.w*(1/4/pi)) + my1.*my2_dy.*my3.*(s.w*(1/4/pi)) + my1.*my2.*my3_dy.*(s.w*(1/4/pi));
my1_dz = -1i*omega*d3./r; my2_dz = 1i*omega*my2.*d3./r; my3_dz = 1./(r.*rr).*(-3*d3./rr.*ddotsn + s.nx(3,:));
Anz_test = my1_dz.*my2.*my3.*(s.w*(1/4/pi)) + my1.*my2_dz.*my3.*(s.w*(1/4/pi)) + my1.*my2.*my3_dz.*(s.w*(1/4/pi));
diff = abs(An - (t.nx(1,:)'.*Anx_test + t.nx(2,:)'.*Any_test + t.nx(3,:)'.*Anz_test));
max(diff(:)) % naive kernel evaluation

%%% redo Anx, separate real and imaginary part
Anx_test2 = (-1i*omega*d1./r).*     exp(1i*omega*r).*                       ddotsn./(r.*rr).*                                   (s.w*(1/4/pi)) ...
          + (1-1i*omega*r).*        (1i*omega*exp(1i*omega*r).*d1./r).*     ddotsn./(r.*rr).*                                   (s.w*(1/4/pi)) ...
          + (1-1i*omega*r).*        exp(1i*omega*r).*                       (1./(r.*rr).*(-3*d1./rr.*ddotsn + s.nx(1,:))).*     (s.w*(1/4/pi));
% Anx_test2 = (1-1i*omega*r).*(1i*omega*exp(1i*omega*r).*d1./r).*ddotsn./(r.*rr).*(s.w*(1/4/pi));
% Anx_test2 = (1-1i*omega*r).*exp(1i*omega*r).*(1./(r.*rr).*(-3*d1./rr.*ddotsn + s.nx(1,:))).*(s.w*(1/4/pi));
% on a fine grid, check smoothness on imaginary part...
d1f = bsxfun(@minus,t.x(1,:)',sf.x(1,:)); % 3 coords of displacement matrix (M*N)
d2f = bsxfun(@minus,t.x(2,:)',sf.x(2,:));
d3f = bsxfun(@minus,t.x(3,:)',sf.x(3,:));
rrf = d1f.^2+d2f.^2+d3f.^2;   % dist^2 mat
rf  = sqrt(rrf);
ddotsnf = bsxfun(@times,d1f,sf.nx(1,:))+bsxfun(@times,d2f,sf.nx(2,:))+bsxfun(@times,d3f,sf.nx(3,:));
Afnx_test2 = (-1i*omega*d1f./rf).*     exp(1i*omega*rf).*                         ddotsnf./(rf.*rrf).*                                      (sf.w*(1/4/pi)) ...
           + (1-1i*omega*rf).*         (1i*omega*exp(1i*omega*rf).*d1f./rf).*     ddotsnf./(rf.*rrf).*                                      (sf.w*(1/4/pi)) ...
           + (1-1i*omega*rf).*         exp(1i*omega*rf).*                         (1./(rf.*rrf).*(-3*d1f./rrf.*ddotsnf + sf.nx(1,:))).*     (sf.w*(1/4/pi));
ui = imag(Anx_test2)*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uif = imag(Afnx_test2)*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
erri = abs(ui-uif)/max(abs(uif(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(erri),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(erri) + " ")

%%% redo Anx, kernel split version, i.e. singularity of lap type times smooth density
lapdlp = ddotsn./(r.*rr);
lapdlpnx = (1./(r.*rr).*(-3*d1./rr.*ddotsn + s.nx(1,:)));
Anx_test3 =  ( lapdlp.*(omega*sin(omega*r)./r.*d1) - 1i*lapdlp.*(omega*cos(omega*r)./r.*d1)) .* (s.w*(1/4/pi)) ... 
           + ( lapdlp.*((-omega*sin(omega*r)+omega^2*r.*cos(omega*r))./r.*d1) + 1i*lapdlp.*( (omega^2*r.*sin(omega*r)+omega*cos(omega*r))./r.*d1 ) ) .* (s.w*(1/4/pi)) ...
           + ( lapdlpnx.*(cos(omega*r)+omega*r.*sin(omega*r)) + 1i*lapdlpnx.*(sin(omega*r)-omega*r.*cos(omega*r)) ) .* (s.w*(1/4/pi)); % check
lapdlpf = ddotsnf./(rf.*rrf);
lapdlpnxf = (1./(rf.*rrf).*(-3*d1f./rrf.*ddotsnf + sf.nx(1,:)));
Afnx_test3 =  ( lapdlpf.*(omega*sin(omega*rf)./rf.*d1f) - 1i*lapdlpf.*(omega*cos(omega*rf)./rf.*d1f)) .* (sf.w*(1/4/pi)) ... 
            + ( lapdlpf.*((-omega*sin(omega*rf)+omega^2*rf.*cos(omega*rf))./rf.*d1f) + 1i*lapdlpf.*( (omega^2*rf.*sin(omega*rf)+omega*cos(omega*rf))./rf.*d1f ) ) .* (sf.w*(1/4/pi)) ...
            + ( lapdlpnxf.*(cos(omega*rf)+omega*rf.*sin(omega*rf)) + 1i*lapdlpnxf.*(sin(omega*rf)-omega*rf.*cos(omega*rf)) ) .* (sf.w*(1/4/pi)); % check
iAnx_test3 =  ( - 1i*lapdlp.*(omega*cos(omega*r)./r.*d1)) .* (s.w*(1/4/pi)) ... 
            + ( 1i*lapdlp.*( (omega^2*r.*sin(omega*r)+omega*cos(omega*r))./r.*d1 ) ) .* (s.w*(1/4/pi)) ...
            + ( 1i*lapdlpnx.*(sin(omega*r)-omega*r.*cos(omega*r)) ) .* (s.w*(1/4/pi)); % check
iAfnx_test3 =  ( - 1i*lapdlpf.*(omega*cos(omega*rf)./rf.*d1f)) .* (sf.w*(1/4/pi)) ... 
            + ( 1i*lapdlpf.*( (omega^2*rf.*sin(omega*rf)+omega*cos(omega*rf))./rf.*d1f ) ) .* (sf.w*(1/4/pi)) ...
            + ( 1i*lapdlpnxf.*(sin(omega*rf)-omega*rf.*cos(omega*rf)) ) .* (sf.w*(1/4/pi)); % check
ui3 = iAnx_test3*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uif3 = iAfnx_test3*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
erri3 = abs(ui3-uif3)/max(abs(uif3(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(erri3),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(erri3) + " ")

%%% Anx, study smoothness of each splitted kernel
% 1. laplace dlp * smooth coeffs
Astar_coefs1 = (omega*sin(omega*r)./r.*d1) .* (s.w*(1/4/pi));
Afstar_coefs1 = (omega*sin(omega*rf)./rf.*d1f) .* (sf.w*(1/4/pi));
ustar1 = Astar_coefs1*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
ufstar1 = Afstar_coefs1*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
errstar1 = abs(ustar1-ufstar1)/max(abs(ufstar1(:)));
figure(2),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(errstar1),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("1st splitted kernel coefs: max abs err " + max(errstar1) + " ")
% 2. 
Astar_coefs2 = ((-omega*sin(omega*r)+omega^2*r.*cos(omega*r))./r.*d1) .* (s.w*(1/4/pi));
Afstar_coefs2 = ((-omega*sin(omega*rf)+omega^2*rf.*cos(omega*rf))./rf.*d1f) .* (sf.w*(1/4/pi));
ustar2 = Astar_coefs2*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
ufstar2 = Afstar_coefs2*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
errstar2 = abs(ustar2-ufstar2)/max(abs(ufstar2(:)));
figure(3),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(errstar2),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("2nd splitted kernel coefs: max abs err " + max(errstar2) + " ")
% 3. 
Astar_coefs3 = (cos(omega*r)+omega*r.*sin(omega*r)) .* (s.w*(1/4/pi));
Afstar_coefs3 = (cos(omega*rf)+omega*rf.*sin(omega*rf)) .* (sf.w*(1/4/pi));
ustar3 = Astar_coefs3*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
ufstar3 = Afstar_coefs3*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
errstar3 = abs(ustar3-ufstar3)/max(abs(ufstar3(:)));
figure(4),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(errstar3),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("3rd splitted kernel coefs: max abs err " + max(errstar3) + " ")

%... done, now do actual laplace dlp and dlpn matrix
AlapD = Lap3dDLPmat(t,s);
diffD = abs(AlapD - lapdlp.*(s.w*(1/4/pi)));
[~,AlapDnx] = Lap3dDLPnmat(t,s);
diffDnx = abs(AlapDnx - lapdlpnx.*(s.w*(1/4/pi)));

Acirc = 1i*imag(Anx_test3);
ADstar = AlapD.*( omega^2*r.*cos(omega*r)./r.*d1 );
ADnxstar = AlapDnx.*(cos(omega*r)+omega*r.*sin(omega*r));
Anx_test4 = Acirc + ADstar + ADnxstar;

keyboard

%%% 0th check, kernel call, then tak imaginary part
A = Helm3dDLPnmat(t,s,omega);
Af = Helm3dDLPnmat(t,sf,omega);
% check if imaginary part is smooth
ui = imag(A)*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
uif = imag(Af)*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';  
% ui = real(A)*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
% uif = real(Af)*fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))';
erri = abs(ui-uif)/max(abs(uif(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(erri),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("max abs err " + max(erri) + " ")

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