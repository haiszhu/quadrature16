% Helmholtz convergence
%
%

% Paths.
folder=fileparts(mfilename('fullpath'));
addpath(fullfile(folder,'..','build'));
addpath(fullfile(folder,'..','utils'));
addpath(fullfile(folder,'..','utils','torus'));

% Layer and wavenumber.
layer='DLP';
% layer='SLP';
zk=1;
ifS=double(strcmp(layer,'SLP')); 
ifD=double(strcmp(layer,'DLP'));

%
maxNumCompThreads(1);
pvalues=2:13; order=16;
tol=eps; maxit=100;
if strcmp(layer,'SLP')
  maxit=1000; % First-kind solve.
  warning_state=warning('off','MATLAB:gmres:tooSmallTolerance');
  warning_cleanup=onCleanup(@() warning(warning_state));
end

%
theta=46*pi/180;
r1=[4*cos(theta);4*sin(theta);0];
r2=[-4*sin(theta);4*cos(theta);0];
g=@(r) exp(1i*zk*vecnorm(r-r1))./vecnorm(r-r1)-exp(1i*zk*vecnorm(r-r2))./vecnorm(r-r2);
s2=linspace(-pi,pi,41);
s2(end)=[];
ctr=[2*cos(s2);2*sin(s2);zeros(size(s2))];
ex=g(ctr);

%
nres=numel(pvalues); Ns=zeros(1,nres); Es=zeros(1,nres); iterations=zeros(1,nres); residuals=zeros(1,nres);
timing=struct('near_mid_seconds',zeros(1,nres),'far_seconds',zeros(1,nres),'setup_seconds',zeros(1,nres),'assembly_seconds',zeros(1,nres),'solve_seconds',zeros(1,nres),'off_surface_seconds',zeros(1,nres),'total_seconds',zeros(1,nres));
fprintf('run_figure4_w16_generic: standard torus, Helmholtz %s, zk=%g%+gi, GL16 density, physical target screening, quadSD_helm_w16 near/self+mid, GL16 far\n',layer,real(zk),imag(zk));
fprintf('    p      N    relL2(U)  near_mid_s    far_s  setup_s assembly_s  solve_s   eval_s  total_s  GMRES\n');

% Fixed 40 interior targets.
off_target = ctr;
off_nh = size(ctr,2);
off_exact = ex(:);

% Far quadrature tables
far_setup_clock = tic;
[g16,w16] = compute_gauss_legendre(order);
P16 = compute_legendre_values(order-1,g16);
w16_tensor = reshape(w16*w16.',1,[])/(2*pi);
timing.far_precompute_seconds = toc(far_setup_clock);
work = zeros(0,1);

% convergence study
for kk=1:nres

  %
  p = pvalues(kk);
  mp = p;
  np = p;
  run_clock = tic;

  %
  order = 16;
  so = [];
  so.mp = mp; % poloidal, fast index, theta
  so.np = np; % toroidal, phi
  so.a = 2; so.b = 1;
  so.orig = -16;
  o.p = order;
  o.root_bound = 1.2;
  setup_clock = tic;
  [pan,N] = create_panels('torus',so,o);
  [sx,snx,~] = getallnodes(pan);
  timing.setup_seconds(kk) = toc(setup_clock);

  % generate exact solution
  alpha = 46*pi/180;
  r1 = [4*cos(alpha); 4*sin(alpha); 0];
  r2 = [-4*sin(alpha); 4*cos(alpha); 0];
  uex = @(x) exp(1i*zk*vecnorm(x-r1,2,1))./vecnorm(x-r1,2,1) - exp(1i*zk*vecnorm(x-r2,2,1))./vecnorm(x-r2,2,1);

  % generate boundary data
  g = uex(sx);

  % generate the naive evaluation matrix
  tx = sx;
  tnx = snx;
  A = complex(zeros(size(sx,2),size(sx,2)));

  % generate near/self and mid quadrature correction
  for k=1:numel(pan)

    s = pan{k};
    t = s.t([2 1],:); % t1: meridian; t2: azimuth
    [i0,j0] = ind2sub([so.mp,so.np],k);
    [i,j] = ind2sub([so.mp,so.np],1:numel(pan));
    offset = [i-i0;j-j0];
    offset = mod(offset+[so.mp;so.np]/2,[so.mp;so.np])-[so.mp;so.np]/2;
    t0_all = repmat(t,1,numel(pan))+repelem(2*offset,1,order^2);

    % Source geometry
    cols = s.inds;
    n = numel(cols);
    setup_clock = tic;

    % Physical target classification
    m = size(tx,2);
    inside16 = false(1,m);
    outside16 = false(1,m);
    [inside16,outside16] = prepare_quadSD_targets_w16_mex(m,tx-s.center,n,s.x-s.center,s.C_X,s.C_Xi,zeros(2,1),ones(2,1),inside16,outside16,zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(1,3),zeros(3,m),0,0,0);
    idxc = find(inside16);
    far = find(outside16);
    m = numel(idxc);
    timing.setup_seconds(kk) = timing.setup_seconds(kk)+toc(setup_clock);

    % Near/self and mid (DLP): Helm3dDLPmat_w16 returns -1/2 of these weights
    near_mid_clock = tic;
    s.C_X_fit = s.fit_src.C_X;
    t = struct('x',tx(:,idxc),'onsurface_comp',true,'otherside_flag',false(1,m),'t0',t0_all(:,idxc),'h',zeros(m,1));
    if ifS
      A(idxc,cols) =  2*Helm3dSLPmat_w16(t,s,zk);
    end
    if ifD
      A(idxc,cols) = -2*Helm3dDLPmat_w16(t,s,zk);
    end
    timing.near_mid_seconds(kk) = timing.near_mid_seconds(kk)+toc(near_mid_clock);

    % k-th source patch far targets
    % actual compute. Normal and area weight from the fitted tangents that near/self and mid use. s.nx and s.w
    % come from D applied to the positions, whose roundoff sets a ~3e-15 floor as patches shrink.
    far_clock = tic;
    X_u = zeros(3,n); X_v = zeros(3,n);
    for coordinate=1:3
      X_u(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,1)*P16.',1,[]);
      X_v(coordinate,:) = reshape(P16*s.C_Xi(:,:,coordinate,2)*P16.',1,[]);
    end
    weighted_normal = cross(X_v,X_u,1).*w16_tensor;   % n*w/(2*pi) at the nodes
    s.nx = weighted_normal./vecnorm(weighted_normal);
    s.w  = 2*pi*vecnorm(weighted_normal);
    if ifS
      A(far,cols) =  2*Helm3dSLPmat(struct('x',tx(:,far)),s,zk);
    end
    if ifD
      A(far,cols) = -2*Helm3dDLPmat(struct('x',tx(:,far)),s,zk);
    end
    timing.far_seconds(kk) = timing.far_seconds(kk)+toc(far_clock);

  end

  timing.assembly_seconds(kk) = toc(run_clock);

  % solve
  rhs = 2*g(:);
  solve_clock = tic;
  [mu,~,residuals(kk),iter] = gmres(@(v) A*v+ifD*v,rhs,[],tol,maxit);
  iterations(kk) = iter(2);
  timing.solve_seconds(kk) = toc(solve_clock);

  % evaluate at targets
  eval_clock = tic;
  [g,wg] = compute_gauss_legendre(26);
  L = compute_legendre_values(15,g)/P16;
  [t1,t2] = ndgrid(g,g);
  tq = [t2(:).';t1(:).'];
  wq = reshape(wg*wg.',1,[]);
  nq = numel(g)^2;
  xq = zeros(3,nq*numel(pan)); Nq = xq;
  muq = complex(zeros(1,nq*numel(pan)));
  for k=1:numel(pan)
    s = pan{k};
    cols = s.inds;
    rows = (k-1)*nq+(1:nq);
    [xq(:,rows),~,~,~,~,Xphi,Xtheta] = s.chart(tq);
    Nq(:,rows) = cross(Xphi,Xtheta,1)*prod(s.psiz);
    muq(rows) = reshape(L*reshape(mu(cols),16,16)*L.',1,[]);
  end

  sq = struct('x',xq,'nx',Nq./vecnorm(Nq),'w',repmat(wq,1,numel(pan)).*vecnorm(Nq));
  if ifS, off_value =  sum(Helm3dSLPmat(struct('x',off_target),sq,zk).*muq,2); end
  if ifD, off_value = -sum(Helm3dDLPmat(struct('x',off_target),sq,zk).*muq,2); end
  timing.off_surface_seconds(kk) = toc(eval_clock);

  %
  Ns(kk) = N;
  Es(kk) = norm(off_value-off_exact)/norm(off_exact);
  timing.total_seconds(kk) = toc(run_clock);
  fprintf('p=%2d N=%6d relL2=%.3e near_mid=%.3fs far=%.3fs setup=%.3fs assembly=%.3fs solve=%.3fs eval=%.3fs total=%.3fs GMRES=%d\n',p,N,Es(kk),timing.near_mid_seconds(kk),timing.far_seconds(kk),timing.setup_seconds(kk),timing.assembly_seconds(kk),timing.solve_seconds(kk),timing.off_surface_seconds(kk),timing.total_seconds(kk),iterations(kk));

end
results = struct('N',Ns,'error',Es,'iterations',iterations,'relres',residuals,'timing',timing);

% convergence plot
figure(4); clf;
loglog(Ns,Es,'b+-','LineWidth',1.2); hold on;
Eref = Es(3)*(Ns/Ns(3)).^(-order/2);
loglog(Ns,Eref,'k--','LineWidth',1);
grid on;
xlabel('Number of discretization points, N');
ylabel('Relative L2 error');
legend('40 fixed interior targets','Order 16: N^{-8}','Location','southwest');
title(sprintf('Torus Helmholtz %s, k=%g%+gi: GL16 density, GL26 evaluation',layer,real(zk),imag(zk)));
% exportgraphics(gcf,fullfile(folder,sprintf('figure4_w16_generic_%s.png',layer)),'Resolution',200);

keyboard
