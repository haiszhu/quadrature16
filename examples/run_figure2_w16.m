% the code is designed to handle nearly singular evaluation in a way
% similar to on-surface, which makes the evaluation routine somewhat
% complicated....
%

folder = fileparts(mfilename('fullpath'));
addpath(fullfile(folder,'..','matlab'))
addpath(fullfile(folder,'..','build'))
addpath(fullfile(folder,'..','utils'))
addpath(fullfile(folder,'..','utils','torus'))

% mp = 8; np = 8;
% mp = 10; np = 10;
mp = 12; np = 12;
% mp = 16; np = 16;

% mp = 12; np = 18;

% select test geometry
order = 16;
geometry = 'torus';
so = [];
so.mp = mp; % poloidal, fast index, theta
so.np = np; % toroidal, phi
so.a = 2; so.b = 1;
so.orig = -16;
o.p = order;
o.root_bound = 1.2;
[pan,N] = create_panels(geometry,so,o);
[sx,snx,sw] = getallnodes(pan);

% generate exact solution
alpha = 46*pi/180;
r1 = [4*cos(alpha); 4*sin(alpha); 0];
r2 = [-4*sin(alpha); 4*cos(alpha); 0];
uex = @(x) 1./vecnorm(x-r1,2,1) - 1./vecnorm(x-r2,2,1);

% generate boundary data
g = uex(sx);

% generate the naive evaluation matrix
tx = sx;
tnx = snx;
% A = Lap3dDLPmat(struct('x',tx),struct('x',sx,'nx',snx,'w',sw));
AS = zeros(size(sx,2),size(sx,2));
AD = zeros(size(sx,2),size(sx,2));

% generate near/self and mid quadrature correction
[g16,w16] = compute_gauss_legendre(order);
P16 = compute_legendre_values(order-1,g16);
w16_tensor = reshape(w16*w16.',1,[])/(2*pi);
work = zeros(0,1);
ifS = 1; ifD = 1; ifSn = 0; ifDn = 0;
% figure(1),clf,
% plot3(sx(1,:),sx(2,:),sx(3,:),'.'); axis equal, hold on
for k=1:numel(pan)
  
  %
  s = pan{k};
  cols = s.inds;
  n = numel(cols);
  setup_clock = tic;

  % near/self and mid region estimate
  m = size(tx,2);
  inside16 = false(1,m);
  outside16 = false(1,m);
  [inside16,outside16] = prepare_quadSD_targets_w16_mex(m,tx-s.center,n,s.x-s.center,s.C_X,s.C_Xi,zeros(2,1),ones(2,1),inside16,outside16,zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(1,3),zeros(3,m),0,0,0);
  idxc = find(inside16);
  far = find(outside16);
  m = numel(idxc);

  % near/self compute
  s.C_X_fit = s.fit_src.C_X;
  t = struct('x',tx(:,idxc),'onsurface_comp',true,'otherside_flag',false(1,m));
  if ifS
    AS(idxc,cols) =  2*Lap3dSLPmat_w16(t,s);
  end
  if ifD
    AD(idxc,cols) = -2*Lap3dDLPmat_w16(t,s);
  end

  % far compute
  X_u = zeros(3,n); X_v = zeros(3,n);
  for coordinate=1:3
    X_u(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,1)*P16.',1,[]);
    X_v(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,2)*P16.',1,[]);
  end
  weighted_normal = cross(X_v,X_u,1).*w16_tensor;   % n*w/(2*pi) at the nodes
  s.nx = weighted_normal./vecnorm(weighted_normal);
  s.w  = 2*pi*vecnorm(weighted_normal);
  if ifS
    AS(far,cols) =  2*Lap3dSLPmat(struct('x',tx(:,far)),s);
  end
  if ifD
    AD(far,cols) = -2*Lap3dDLPmat(struct('x',tx(:,far)),s);
  end

  % keyboard
end

% gmres solve
constant_error = norm(AD*ones(N,1)-1,inf);
[mu,flag,relres,iter] = gmres(@(v) v+AD*v,2*g(:),[],eps,100);

% Interior target grid
tx_scale = 2;
xgrid = linspace(so.a-tx_scale*so.b,so.a+tx_scale*so.b,21);
ygrid = linspace(-tx_scale*so.b,tx_scale*so.b,21);
zgrid = linspace(-so.b,so.b,21);
[X,Y,Z] = meshgrid(xgrid,ygrid,zgrid);
margin = 1e-4*so.b;
inside = (hypot(X,Y)-so.a).^2+Z.^2 < (so.b-margin)^2;
targs = [X(inside).';Y(inside).';Z(inside).'];

% evaluate solution
u = zeros(size(targs,2),1);
exact = uex(targs).';

% via quadSD
for k=1:numel(pan)
  %
  s = pan{k};
  cols = s.inds;
  n = numel(cols);
  setup_clock = tic;

  % near/self and mid region rough estimate
  m = size(targs,2);
  inside16 = false(1,m);
  outside16 = false(1,m);
  [inside16,outside16] = prepare_quadSD_targets_w16_mex(m,targs-s.center,n,s.x-s.center,s.C_X,s.C_Xi,zeros(2,1),ones(2,1),inside16,outside16,zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(1,3),zeros(3,m),0,0,0);

  %
  idxc = find(inside16);
  far = find(outside16);
  m = numel(idxc);

  % near/self compute
  s.C_X_fit = s.fit_src.C_X;
  t = struct('x',targs(:,idxc),'onsurface_comp',false,'otherside_flag',false(1,m));
  u(idxc) = u(idxc) - Lap3dDLPmat_w16(t,s)*mu(cols);

  % far compute
  X_u = zeros(3,n); X_v = zeros(3,n);
  for coordinate=1:3
    X_u(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,1)*P16.',1,[]);
    X_v(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,2)*P16.',1,[]);
  end
  weighted_normal = cross(X_v,X_u,1).*w16_tensor;   % n*w/(2*pi) at the nodes
  s.nx = weighted_normal./vecnorm(weighted_normal);
  s.w  = 2*pi*vecnorm(weighted_normal);
  u(far) = u(far) - Lap3dDLPmat(struct('x',targs(:,far)),s)*mu(cols);

  % keyboard
end

%
abs_error = abs(u-exact);
[max_error,worst] = max(abs_error);
rel_error = norm(u-exact,inf)/norm(exact,inf);
fprintf('Interior targets: %d; max error %.3e; relative Linf %.3e\n',numel(u),max_error,rel_error);
fprintf('Worst target: (%.16g, %.16g, %.16g)\n',targs(:,worst));
error_grid = nan(size(X));
error_grid(inside) = abs_error;

figure(3),clf
scatter3(targs(1,:),targs(2,:),targs(3,:),18,log10(max(abs_error,1e-16)),'filled');
hold on;
tb = linspace(0,2*pi,201);
for k = 1:numel(pan)
  xb = pan{k}.bdchart(tb);
  plot3(xb(1,:),xb(2,:),xb(3,:),'r-','LineWidth',1);
end
axis equal;
view(3);
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
cb = colorbar;
cb.Label.String = 'log_{10}|u-u_{exact}|';
title('Interior evaluation error');

keyboard