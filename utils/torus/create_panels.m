function [pan N] = create_panels(shape,so,o)
% CREATE_PANELS  Build panels for various 3D surface geometries
%
% [pan N] = create_panels(shape,shapeopts,convopts)
%  Creates struct array of panels with their quadrature nodes for unknowns
%
%  shape = 'torus': shapeopts has fields a (major radius), b (minor rad),
%                   and optional fields:
%                     mp (# panels small circ),
%                     np (# pan big circ),
%                     orig = 0 (panels start at phi=th=0)
%                            1 (1/2-panel offsets, default)
%                   mp * np panels are ordered in incr theta (i) fast,
%                   incr phi (j) slow  (ie flipud of matrix ordering).
%                   Dof ordering within panels is same (y fast, x slow).
%                   If b is a cell array of 3 function handles, interpreted as
%                   { b(theta,phi), b_theta(theta,phi), and b_phi(theta,phi },
%                   a minor radius function over theta and phi in [0,2pi), and
%                   its correct partials. This allows height-modulated torus.
%
% Output
%  N : total number of unknowns
%  pan : cell array with each element having fields
%        (n = p^2 is # nodes per panel):
%        inds - indices from full N list belonging in this panel
%        t - (2 by n) chart parameters of smooth nodes
%        x - (3 by n) 3D coords of smooth nodes
%        nx - (3 by n) 3D coords of normal derivs at smooth nodes
%        w - (1 by n) smooth quad weights including speed factors
%        [sp - (1 by n) speed factors at nodes (Jacobean area factor)]
%        chart - parametrization function for panels, specifically:
%                input 2-by-N pts in [-1,1]^2, output [x n_x speed] where
%                x and nx are 3-by-N (pts and normals in R3), speed is 1-by-N.
%                (note speed is relative to global param of shape)
%        spfac - (scalar) size factor converting to speed in chart coords.
%        nei - indices in panel list of (usually 8) neighbor panels, in a
%              standard ordering.
%        neix1 - (3 by 8) 3D coords of first node in each neighbor panel,
%              in std ordering, useful to identify neighbors w/o global index
%  In addition pan{1} contains the following fields needed to specify geometry:
%        mp, np - numbers of vertical and horizontal panels.
%        topo - topology, eg 'torus'
%
% Optional GL16 torus preparation:
% o.root_bound - extended sampling limit
% center, C_X, C_Xi, C_Xij, fit_src
% x_ext, t_ext, C_Xi_ext, C_Xij_ext, fit_ext, root_bound

% Barnett 7/11/15. restarted 7/14/16, swapped np & mp 7/18/16
% 12/29/16 neix1. torus allowed wobbly 10/10/18, so.orig 12/12/18

% future ideas:
% *** why not make just the surf pts be input, and derive all geom
% from that? But, let the q used for that differ from the p used for
% density rep.
%   *** Use bkwd stable monomial solve to get the local charts for
%  r, r_u, r_v, n, etc ?

% Hai 08/22/22, panel_smooth_quad compute nx & sp numerically from
% pan{k}.x; add panelbd_smooth_quad for panel boundary; add pan{k}.bdchart
% and pan{k}.ltrichart & pan{k}.ltrichart

if nargin==0, test_create_panels; return; end
if nargin<3, o=[]; end
if isfield(o,'p'), p = o.p; else p = 16; end  % # nodes per side of panel (order)
if isfield(o,'fquad'), fquad = o.fquad; else fquad = 1; end 

rboxZ = @(t) (mod(t,2*pi)<pi/2).*(-1+4/pi*mod(t,2*pi)-1i)...
            +((mod(t,2*pi)>=pi/2)&(mod(t,2*pi)<pi)).*(1+1i*(-1+4/pi*(mod(t,2*pi)-pi/2)))...
            +((mod(t,2*pi)>=pi)&(mod(t,2*pi)<3*pi/2)).*(1-4/pi*(mod(t,2*pi)-pi)+1i)...
            +(mod(t,2*pi)>=3*pi/2).*(-1+1i*(1-4/pi*(mod(t,2*pi)-3*pi/2)));
striZ = @(t) (mod(t,2*pi)<pi/2).*(0+1i*(1-2/pi*mod(t,2*pi)))...
            +((mod(t,2*pi)>=pi/2)&(mod(t,2*pi)<pi)).*(2/pi*(mod(t,2*pi)-pi/2)+0*1i)...
            +(mod(t,2*pi)>=pi).*(1-1/pi*(mod(t,2*pi)-pi)+1i*(1/pi*(mod(t,2*pi)-pi)));

if strcmp(shape,'torus') || strcmp(shape,'chestnut')  || strcmp(shape,'stellarator') 
  if isfield(so,'mp'), mp = so.mp; else mp = 8; end  % # pans on small circle
  if isfield(so,'np'), np = so.np; else np = 12; end  % # pans on big circle
  if isfield(so,'orig'), orig = so.orig; 
  else
    if strcmp(shape,'torus'), orig = 1; end
    if strcmp(shape,'chestnut'), orig = 0; end
    if strcmp(shape,'stellarator'), orig = 1; end
  end  % pan origin shift
  tpansiz1 = pi/np; tpansiz2 = pi/mp;     % panel half-sizes in parameter
  npan = mp*np; N = npan*p^2;       % same order all panels, total # nodes
  panind = reshape(1:npan,[mp np]); % 2d grid of panel numbers (sets dof order)
  [xg,~] = compute_gauss_legendre(p); [xg1,xg2] = meshgrid(xg); xx = [xg1(:)';xg2(:)']; step = 1e-30; % GL nodes, as in panel_smooth_quad
  pan = cell(npan,1);
  for j=1:np, for i=1:mp, k=panind(i,j);     % loop over panels (in any order)
      pan{k}.psiz = [tpansiz1,tpansiz2];     % param sizes
      pan{k}.inds = p^2*(k-1)+(1:p^2);     % set up indices of unknowns
      switch shape
        case 'torus'
          pan{k}.chart = @(t) torusparam(so.a,so.b,(t(1,:)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2);   % note speed is relative to original params; note i up, j across
          pan{k}.hcurvature = @(t) toruscurvature(so.a,so.b,(t(1,:)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2); 
          % Nodal geometry: first derivatives from the chart, second by complex step.
          [~,~,~,~,~,Xphi,Xtheta] = pan{k}.chart(xx);
          tu = xx; tu(2,:) = tu(2,:)+1i*step;
          tv = xx; tv(1,:) = tv(1,:)+1i*step;
          [~,~,~,~,~,Xphiu,Xthetau] = pan{k}.chart(tu);
          [~,~,~,~,~,Xphiv,~] = pan{k}.chart(tv);
          pan{k}.Xu = tpansiz2*Xtheta;
          pan{k}.Xv = tpansiz1*Xphi;
          pan{k}.Xuu = tpansiz2*imag(Xthetau)/step;
          pan{k}.Xuv = tpansiz1*imag(Xphiu)/step;
          pan{k}.Xvv = tpansiz1*imag(Xphiv)/step;
        case 'chestnut'
          pan{k}.chart = @(t) chestnutparam(so.a,so.b,(t(1,:)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2);
        case 'stellarator'
          pan{k}.chart = @(t) stellaratorparam((t(1,:)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2); 
      end
      pan{k}.spfac = tpansiz1*tpansiz2;  % factor to convert speed to chart z
      nei = panind(mod([i-2 i-1 i],mp)+1,mod([j-2 j-1 j],np)+1); % ordering key
      nei = nei(nei~=k);    % remove self from list, to leave 8
      pan{k}.nei = nei(:);
      % extra: bd + tri
      switch shape
        case 'torus'
          pan{k}.bdchart = @(t) torusparam(so.a,so.b,(real(rboxZ(t))+orig-1+2*j)*tpansiz1,(imag(rboxZ(t))+orig-1+2*i)*tpansiz2);  % t in [0,2*pi]
          pan{k}.ltrichart = @(uvs) torusparam(so.a,so.b,( 2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,( 2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.utrichart = @(uvs) torusparam(so.a,so.b,(-2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,(-2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.ltriuvs = @(t) [( t(1,:)/tpansiz1-orig+1-2*j)/2+1/2;  (t(2,:)/tpansiz2-orig+1-2*i)/2+1/2]; % from parametric space [0,2*pi]x[0,2*pi] to chart uvs
          pan{k}.utriuvs = @(t) [-(t(1,:)/tpansiz1-orig+1-2*j)/2+1/2; -(t(2,:)/tpansiz2-orig+1-2*i)/2+1/2];
          pan{k}.bdltrichart = @(t) torusparam(so.a,so.b,( 2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,( 2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
          pan{k}.bdutrichart = @(t) torusparam(so.a,so.b,(-2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,(-2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
        case 'chestnut'
          pan{k}.bdchart = @(t) chestnutparam(so.a,so.b,(real(rboxZ(t))+orig-1+2*j)*tpansiz1,(imag(rboxZ(t))+orig-1+2*i)*tpansiz2);  % t in [0,2*pi]
          pan{k}.ltrichart = @(uvs) chestnutparam(so.a,so.b,( 2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,( 2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.utrichart = @(uvs) chestnutparam(so.a,so.b,(-2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,(-2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.ltriuvs = @(t) [( t(1,:)/tpansiz1-orig+1-2*j)/2+1/2;  (t(2,:)/tpansiz2-orig+1-2*i)/2+1/2]; % from parametric space [0,2*pi]x[0,2*pi] to chart uvs
          pan{k}.utriuvs = @(t) [-(t(1,:)/tpansiz1-orig+1-2*j)/2+1/2; -(t(2,:)/tpansiz2-orig+1-2*i)/2+1/2];
          pan{k}.bdltrichart = @(t) chestnutparam(so.a,so.b,( 2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,( 2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
          pan{k}.bdutrichart = @(t) chestnutparam(so.a,so.b,(-2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,(-2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
        case 'stellarator'
          pan{k}.bdchart = @(t) stellaratorparam((real(rboxZ(t))+orig-1+2*j)*tpansiz1,(imag(rboxZ(t))+orig-1+2*i)*tpansiz2);  % t in [0,2*pi]
          pan{k}.ltrichart = @(uvs) stellaratorparam(( 2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,( 2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.utrichart = @(uvs) stellaratorparam((-2*(uvs(1,:)-1/2)+orig-1+2*j)*tpansiz1,(-2*(uvs(2,:)-1/2)+orig-1+2*i)*tpansiz2); % uvs from get_vioreanu_nodes(p)
          pan{k}.ltriuvs = @(t) [( t(1,:)/tpansiz1-orig+1-2*j)/2+1/2;  (t(2,:)/tpansiz2-orig+1-2*i)/2+1/2]; % from parametric space [0,2*pi]x[0,2*pi] to chart uvs
          pan{k}.utriuvs = @(t) [-(t(1,:)/tpansiz1-orig+1-2*j)/2+1/2; -(t(2,:)/tpansiz2-orig+1-2*i)/2+1/2];
          pan{k}.bdltrichart = @(t) stellaratorparam(( 2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,( 2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
          pan{k}.bdutrichart = @(t) stellaratorparam((-2*(real(striZ(t))-1/2)+orig-1+2*j)*tpansiz1,(-2*(imag(striZ(t))-1/2)+orig-1+2*i)*tpansiz2);
      end
      % extra: smooth quadr (debug)
      if fquad, pan{k}.panel_smooth_quadk = @(p) panel_smooth_quadk(pan{k}.chart,pan{k}.psiz,p); end
      % extra: panel upper and lower limits
      pan{k}.tlo = [(-1+orig-1+2*j)*tpansiz1; (-1+orig-1+2*i)*tpansiz2];
      pan{k}.thi = [( 1+orig-1+2*j)*tpansiz1; ( 1+orig-1+2*i)*tpansiz2];
    end,end
    pan{1}.mp=mp; pan{1}.np=np; pan{1}.topo = 'torus'; pan{1}.p = p;
end
pan = panel_smooth_quad(pan,p);  % does all panels
pan = panelbd_smooth_quad(pan,p); % does all panels' boundaries
if isfield(o,'root_bound')
  assert(strcmp(shape,'torus') && p==16,'create_panels:preparation','Geometry preparation requires GL16 torus panels.');
  assert(isnumeric(o.root_bound) && isreal(o.root_bound) && isscalar(o.root_bound) && isfinite(o.root_bound) && o.root_bound>=1,'create_panels:rootBound','root_bound must be a finite real scalar >= 1.');
  g_ext = o.root_bound*xg;
  for k=1:numel(pan)
    [pan{k}.center,pan{k}.C_X,pan{k}.C_Xi,pan{k}.C_Xij,pan{k}.fit_src] = cxinit_3d(xg,pan{k}.x,pan{k}.Xu,pan{k}.Xv,pan{k}.Xuu,pan{k}.Xuv,pan{k}.Xvv);
    pan{k}.root_bound = o.root_bound;
    pan{k}.t_ext = o.root_bound*pan{k}.t;
    [pan{k}.x_ext,~,~,~,~,Xphi,Xtheta] = pan{k}.chart(pan{k}.t_ext);
    tu = pan{k}.t_ext; tu(2,:) = tu(2,:)+1i*step;
    tv = pan{k}.t_ext; tv(1,:) = tv(1,:)+1i*step;
    [~,~,~,~,~,Xphiu,Xthetau] = pan{k}.chart(tu);
    [~,~,~,~,~,Xphiv,~] = pan{k}.chart(tv);
    Xu_ext = pan{k}.psiz(2)*Xtheta;
    Xv_ext = pan{k}.psiz(1)*Xphi;
    Xuu_ext = pan{k}.psiz(2)*imag(Xthetau)/step;
    Xuv_ext = pan{k}.psiz(1)*imag(Xphiu)/step;
    Xvv_ext = pan{k}.psiz(1)*imag(Xphiv)/step;
    [~,~,pan{k}.C_Xi_ext,pan{k}.C_Xij_ext,pan{k}.fit_ext] = cxinit_3d(g_ext,pan{k}.x_ext,Xu_ext,Xv_ext,Xuu_ext,Xuv_ext,Xvv_ext);
  end
end
% add neix1 tags...
for k=1:numel(pan), nei=pan{k}.nei;
  x = []; for n=1:numel(nei), x = [x pan{nei(n)}.x(:,1)]; end
  pan{k}.neix1 = x;
end
%%%

function pan = panelbd_smooth_quad(pan,p)   % setup G-L boundary native quadr
np = 4; x = zeros(p*np,1); w = x;
[x0, w0, D] = compute_gauss_legendre(p); 
xlo = (0:np-1)/np*2-1; % panel start params
xhi = (1:np)/np*2-1; % panel end params
pt = xhi - xlo;   
for i=1:np
    ii = (i-1)*p+(1:p); % indices of this panel
    x(ii) = xlo(i) + (1+x0')/2*pt(i); 
    w(ii) = w0*pt(i)/2;
end
x1 = [x;ones(size(x));-x;-ones(size(x));-x;x]; 
x2 = [-ones(size(x));x;ones(size(x));-x;x;-x];
xx = [x1(:)';x2(:)']; ww = w(:)';
xp = zeros(3,np*p);
for m=1:numel(pan)
  xbd = pan{m}.chart(xx);  % get 1 output
  for k=1:6     % 6 boundaries (two triangles)
    for i=1:np
      ii = (i-1)*p+(1:p); % indices of this panel
      xp(1,ii) = D*xbd(1,(k-1)*p*np+ii)'*2/pt(i); xp(2,ii) = D*xbd(2,(k-1)*p*np+ii)'*2/pt(i); xp(3,ii) = D*xbd(3,(k-1)*p*np+ii)'*2/pt(i);
    end
    sp = sqrt(xp(1,:).^2+xp(2,:).^2+xp(3,:).^2);
    pan{m}.tang{k} = xp./(ones(3,1)*sp);
    pan{m}.xbd{k} = xbd(:,(k-1)*p*np+(1:p*np));
    pan{m}.ws{k} = ww.*sp;
  end
end


function pan = panel_smooth_quad(pan,p)   % setup G-L tensor native quadr
[x, w, D] = compute_gauss_legendre(p);   % Gauss-Legendre nodes, weights, differentiation matrix on [-1,1]
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
ww = w(:)*w(:)'; ww = ww(:)';                     % 1*p^2 weights
for m=1:numel(pan)
  % only use pan{m}.x
  if 0
    [pan{m}.x pan{m}.nx sp] = pan{m}.chart(xx);  % get 1
  else
    pan{m}.x = pan{m}.chart(xx);  % get 3
    % compute pan{m}.nx & sp from pan{m}.x
    r = pan{m}.x;
    tp_ind = reshape(1:p^2,[p p]);
    xts = nan(1,p^2); yts = xts; zts = xts; % partial t (small circle)
    for k=1:p
      xk = r(1,tp_ind(:,k)); yk = r(2,tp_ind(:,k)); zk = r(3,tp_ind(:,k));
      xts(tp_ind(:,k)) = 1/pan{m}.psiz(2)*D*xk(:); yts(tp_ind(:,k)) = 1/pan{m}.psiz(2)*D*yk(:); zts(tp_ind(:,k)) = 1/pan{m}.psiz(2)*D*zk(:);
    end
    rts = [xts(:),yts(:),zts(:)]';
    xps = nan(1,p^2); yps = xps; zps = xps; % partial p (big circle)
    for k=1:p
      x = r(1,tp_ind(k,:)); y = r(2,tp_ind(k,:)); z = r(3,tp_ind(k,:));
      xps(tp_ind(k,:)) = 1/pan{m}.psiz(1)*D*x(:); yps(tp_ind(k,:)) = 1/pan{m}.psiz(1)*D*y(:); zps(tp_ind(k,:)) = 1/pan{m}.psiz(1)*D*z(:);
    end
    rps = [xps(:),yps(:),zps(:)]';
    dt = sqrt(sum(rts.^2,1));
    dp = sqrt(sum(rps.^2,1));    % could reuse this data
    nx =  cross(rps, rts);       % outward normal
    sp = sqrt(sum(nx.^2,1)); % speeds
    pan{m}.nx = nx.*repmat(1./sp, [3 1]);
  end

  pan{m}.w = sp.*ww*pan{m}.spfac;
  pan{m}.N = p*p;
  pan{m}.t = xx;          % same for all panels, for now
  
end

%%%%%%%%%%%%%%%%%%%%%%
function test_create_panels  % self-test
addpath(fileparts(fileparts(mfilename('fullpath'))))
type = 'torus', so.a = 1; so.b = 0.5; o = [];
[pan N] = create_panels(type,so,o);

% plot it...
figure; x = getallnodes(pan);
plot3(x(1,:),x(2,:),x(3,:),'.','markersize',1); axis equal vis3d
k = 57;    % random panel (conveniently faces viewer for np=8, mp=12)
hold on; for i=1:8, knei=pan{k}.nei(i);
  x = pan{knei}.x; plot3(x(1,:),x(2,:),x(3,:),'g.','markersize',5);
  l = mean(x,2); text(l(1),l(2),l(3),sprintf('%d(%d)',knei,i+1),'fontsize',20)
end
hold on; x = pan{k}.x; plot3(x(1,:),x(2,:),x(3,:),'r.','markersize',10); % self
l = mean(x,2); text(l(1),l(2),l(3),sprintf('%d(%d)',k,1),'fontsize',20,'color',[1 0 0])
l = 0.1; nx = pan{k}.nx; plot3([x(1,:);x(1,:)+l*nx(1,:)],[x(2,:);x(2,:)+l*nx(2,:)],[x(3,:);x(3,:)+l*nx(3,:)],'k-');  % normals
text(x(1,1),x(2,1),x(3,1),'dof 1'); text(x(1,2),x(2,2),x(3,2),'dof 2');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
pan{k}   % look inside struct
%pan{k}.nei

disp('surface area test:')
np = 3; so.np = np; so.mp = ceil(1.5*np); 
[pan N] = create_panels(type,so,o);
[~,~,w] = getallnodes(pan);    % get all weights
fprintf('np=%d: area err = %.3g\n',np,sum(w) - 2*pi^2)

disp('Gauss law flux test:')
zo = [0.9; -0.2; 0.1];    % src pt, must be inside the shape
plot3(zo(1),zo(2),zo(3),'k.','markersize',20);
for mp = 4:2:12, so.mp = mp; so.np = ceil(1.5*mp); 
  [pan N] = create_panels(type,so,o);
  [x,nx,w] = getallnodes(pan);
  d = bsxfun(@minus,x,zo); r = sqrt(sum(d.^2,1));
  ddotn = sum(d.*nx,1);
  fprintf('mp=%d (N=%d): flux err = %.3g\n',mp,N,sum(w.*ddotn./r.^3)/(4*pi) - 1.0)
end
