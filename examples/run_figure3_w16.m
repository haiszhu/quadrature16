% the code is designed to handle nearly singular evaluation in a way
% similar to on-surface, which makes the evaluation routine somewhat
% complicated....
%
% Stokes

folder = fileparts(mfilename('fullpath'));
addpath(fullfile(folder,'..','matlab'))
addpath(fullfile(folder,'..','build'))
addpath(fullfile(folder,'..','utils'))
addpath(fullfile(folder,'..','utils','torus'))

% mp = 4; np = 4;
mp = 6; np = 6;
% mp = 8; np = 8;
mp = 10; np = 10;

% select test geometry
order = 16;
so = [];
so.mp = mp; % poloidal, fast index, theta
so.np = np; % toroidal, phi
so.a = 2; so.b = 1;
so.orig = -16;
o.p = order;
o.root_bound = 1.2;
[pan,N] = create_panels('torus',so,o);
[sx,snx,sw] = getallnodes(pan);

% generate exact solution: two Stokeslets outside the torus (run_figure3_w16_generic.m)
alpha = 46*pi/180;
r1 = [4*cos(alpha); 4*sin(alpha); 0];
r2 = [-4*sin(alpha); 4*cos(alpha); 0];
forces = [1,-.4;.3,.7;-.2,.5];
uex = @(x) reshape(Sto3dSLPmat(struct('x',x),struct('x',[r1,r2],'w',[1,1]))*reshape(forces.',[],1),size(x,2),3).';

% generate boundary data
g = uex(sx);   % 3 x N

% D in quadSD_stok_w16's convention: -Sto3dDLPmat_w16 near/self and mid, -Sto3dDLPmat far
A = zeros(3*N,3*N);
[g16,w16] = compute_gauss_legendre(order);
P16 = compute_legendre_values(order-1,g16);
w16_tensor = reshape(w16*w16.',1,[])/(2*pi);
for k=1:numel(pan)

  %
  s = pan{k};
  cols = s.inds;
  n = numel(cols);
  cols3 = [cols,cols+N,cols+2*N];

  % near/self and mid region estimate
  m = size(sx,2);
  inside16 = false(1,m);
  outside16 = false(1,m);
  [inside16,outside16] = prepare_quadSD_targets_w16_mex(m,sx-s.center,n,s.x-s.center,s.C_X,s.C_Xi,zeros(2,1),ones(2,1),inside16,outside16,zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(1,3),zeros(3,m),0,0,0);
  idxc = find(inside16);
  far = find(outside16);
  m = numel(idxc);

  % near/self compute
  s.C_X_fit = s.fit_src.C_X;
  t = struct('x',sx(:,idxc),'onsurface_comp',true,'otherside_flag',false(1,m));
  A([idxc,idxc+N,idxc+2*N],cols3) = -Sto3dDLPmat_w16(t,s);

  % far compute, with the normal and area weight of the fitted tangents that near/self and mid use
  X_u = zeros(3,n); X_v = zeros(3,n);
  for coordinate=1:3
    X_u(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,1)*P16.',1,[]);
    X_v(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,2)*P16.',1,[]);
  end
  weighted_normal = cross(X_v,X_u,1).*w16_tensor;   % n*w/(2*pi) at the nodes
  s.nx = weighted_normal./vecnorm(weighted_normal);
  s.w  = 2*pi*vecnorm(weighted_normal);
  A([far,far+N,far+2*N],cols3) = -Sto3dDLPmat(struct('x',sx(:,far)),s);

  % keyboard
end

% gmres solve: 1/2*mu+D*mu = g with the rank-one completion with the normal of run_figure3_w16_generic.m
normal = reshape(snx.',[],1);
normal_w = reshape((snx.*sw).',[],1);
area = sum(sw);
applyA = @(v) A*v+normal*(normal_w.'*v)/area;
rhs = reshape(g.',[],1);
[mu,iter,~,relres] = myGMRES(applyA,rhs,3*N,100,eps,0.5);

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
u = zeros(3,size(targs,2));
exact = uex(targs);

% via quadSD
for k=1:numel(pan)
  %
  s = pan{k};
  cols = s.inds;
  n = numel(cols);
  cols3 = [cols,cols+N,cols+2*N];

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
  u(:,idxc) = u(:,idxc) - reshape(Sto3dDLPmat_w16(t,s)*mu(cols3),m,3).';

  % far compute, same fitted normal and weight as the solve
  X_u = zeros(3,n); X_v = zeros(3,n);
  for coordinate=1:3
    X_u(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,1)*P16.',1,[]);
    X_v(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,2)*P16.',1,[]);
  end
  weighted_normal = cross(X_v,X_u,1).*w16_tensor;
  s.nx = weighted_normal./vecnorm(weighted_normal);
  s.w  = 2*pi*vecnorm(weighted_normal);
  u(:,far) = u(:,far) - reshape(Sto3dDLPmat(struct('x',targs(:,far)),s)*mu(cols3),numel(far),3).';

  % keyboard
end

%
abs_error = vecnorm(u-exact,2,1);
[max_error,worst] = max(abs_error);
rel_error = max_error/max(vecnorm(exact,2,1));
fprintf('Stokes DLP, %dx%d patches: GMRES %d iterations, relres %.1e\n',mp,np,iter,relres);
fprintf('Interior targets: %d; max error %.3e; relative Linf %.3e\n',numel(abs_error),max_error,rel_error);
fprintf('Worst target: (%.16g, %.16g, %.16g)\n',targs(:,worst));

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
title('Stokes interior evaluation error');

keyboard
