% Stokes convergence (stokes full matrix requires larger ram)
%
%

% Paths.
folder=fileparts(mfilename('fullpath'));
addpath(fullfile(folder,'..','build'));
addpath(fullfile(folder,'..','utils'));
addpath(fullfile(folder,'..','utils','torus'));

% Layer.
layer='DLP';
% layer='SLP';
ifS=double(strcmp(layer,'SLP')); 
ifD=double(strcmp(layer,'DLP'));

%
maxNumCompThreads(1);
pvalues=2:8; order=16;
tol=eps; maxit=100;
if ifS, maxit=1000; end % First-kind solve.

% Exact solution: two Stokeslets outside the torus.
theta=46*pi/180;
r1=[4*cos(theta);4*sin(theta);0];
r2=[-4*sin(theta);4*cos(theta);0];
forces=[1,-.4;.3,.7;-.2,.5];
uex=@(x) reshape(Sto3dSLPmat(struct('x',x),struct('x',[r1,r2],'w',[1,1]))*reshape(forces.',[],1),size(x,2),3).';
s2=linspace(-pi,pi,41);
s2(end)=[];
ctr=[2*cos(s2);2*sin(s2);zeros(size(s2))];

%
nres=numel(pvalues); Ns=zeros(1,nres); Es=zeros(1,nres); iterations=zeros(1,nres); residuals=zeros(1,nres);
timing=struct('near_mid_seconds',zeros(1,nres),'far_seconds',zeros(1,nres),'setup_seconds',zeros(1,nres),'assembly_seconds',zeros(1,nres),'solve_seconds',zeros(1,nres),'off_surface_seconds',zeros(1,nres),'total_seconds',zeros(1,nres));
fprintf('run_figure3_w16_generic: standard torus, Stokes %s, GL16 density, physical target screening, quadSD_stok_w16 near/self+mid, GL16 far\n',layer);

% Fixed 40 interior targets; velocity stacked [u1;u2;u3].
off_target = ctr;
off_exact = reshape(uex(ctr).',[],1);

% Far quadrature tables
far_setup_clock = tic;
[g16,w16] = compute_gauss_legendre(order);
P16 = compute_legendre_values(order-1,g16);
w16_tensor = reshape(w16*w16.',1,[])/(2*pi);
timing.far_precompute_seconds = toc(far_setup_clock);

% convergence study
for kk=1:nres

  %
  p = pvalues(kk);
  mp = p;
  np = p;
  run_clock = tic;

  %
  so = [];
  so.mp = mp; % poloidal, fast index, theta
  so.np = np; % toroidal, phi
  so.a = 2; so.b = 1;
  so.orig = -16;
  o.p = order;
  o.root_bound = 1.2;
  setup_clock = tic;
  [pan,N] = create_panels('torus',so,o);
  [sx,snx,sw] = getallnodes(pan);
  timing.setup_seconds(kk) = toc(setup_clock);

  % generate boundary data, 3 x N
  g = uex(sx);

  % generate the naive evaluation matrix; unknowns stacked [mu1;mu2;mu3]
  tx = sx;
  A = zeros(3*N,3*N);

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
    cols3 = [cols,cols+N,cols+2*N];
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

    % Near/self and mid: SLP as Sto3dSLPmat; DLP in quadSD_stok_w16's convention (r = y-x), the negative of Sto3dDLPmat_w16
    near_mid_clock = tic;
    s.C_X_fit = s.fit_src.C_X;
    t = struct('x',tx(:,idxc),'onsurface_comp',true,'otherside_flag',false(1,m),'t0',t0_all(:,idxc),'h',zeros(m,1));
    if ifS
      A([idxc,idxc+N,idxc+2*N],cols3) =  Sto3dSLPmat_w16(t,s);
    end
    if ifD
      A([idxc,idxc+N,idxc+2*N],cols3) = -Sto3dDLPmat_w16(t,s);
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
      A([far,far+N,far+2*N],cols3) =  Sto3dSLPmat(struct('x',tx(:,far)),s);
    end
    if ifD
      A([far,far+N,far+2*N],cols3) = -Sto3dDLPmat(struct('x',tx(:,far)),s);
    end
    timing.far_seconds(kk) = timing.far_seconds(kk)+toc(far_clock);

  end

  timing.assembly_seconds(kk) = toc(run_clock);

  % solve: SLP S*mu = g (S[n] = 0 on a closed surface), DLP 1/2*mu+D*mu = g; both with a rank-one
  % completion with the normal.
  normal = reshape(snx.',[],1);
  normal_w = reshape((snx.*sw).',[],1);
  area = sum(sw);
  applyA = @(v) A*v+normal*(normal_w.'*v)/area;
  rhs = reshape(g.',[],1);
  solve_clock = tic;
  [mu,iterations(kk),~,residuals(kk)] = myGMRES(applyA,rhs,3*N,maxit,tol,0.5*ifD);
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
  muq = zeros(3,nq*numel(pan));
  for k=1:numel(pan)
    s = pan{k};
    cols = s.inds;
    rows = (k-1)*nq+(1:nq);
    [xq(:,rows),~,~,~,~,Xphi,Xtheta] = s.chart(tq);
    Nq(:,rows) = cross(Xphi,Xtheta,1)*prod(s.psiz);
    for ic=1:3
      muq(ic,rows) = reshape(L*reshape(mu(cols+(ic-1)*N),16,16)*L.',1,[]);
    end
  end

  sq = struct('x',xq,'nx',Nq./vecnorm(Nq),'w',repmat(wq,1,numel(pan)).*vecnorm(Nq));
  if ifS, off_value =  sum(Sto3dSLPmat(struct('x',off_target),sq).*reshape(muq.',1,[]),2); end
  if ifD, off_value = -sum(Sto3dDLPmat(struct('x',off_target),sq).*reshape(muq.',1,[]),2); end
  timing.off_surface_seconds(kk) = toc(eval_clock);

  %
  Ns(kk) = N;
  Es(kk) = norm(off_value-off_exact)/norm(off_exact);
  timing.total_seconds(kk) = toc(run_clock);
  fprintf('p=%2d N=%6d relL2=%.3e near_mid=%.3fs far=%.3fs setup=%.3fs assembly=%.3fs solve=%.3fs eval=%.3fs total=%.3fs GMRES=%d\n',p,N,Es(kk),timing.near_mid_seconds(kk),timing.far_seconds(kk),timing.setup_seconds(kk),timing.assembly_seconds(kk),timing.solve_seconds(kk),timing.off_surface_seconds(kk),timing.total_seconds(kk),iterations(kk));

end
results = struct('N',Ns,'error',Es,'iterations',iterations,'relres',residuals,'timing',timing);

% convergence plot
figure(1); clf;
loglog(Ns,Es,'b+-','LineWidth',1.2); hold on;
Eref = Es(3)*(Ns/Ns(3)).^(-order/2);
loglog(Ns,Eref,'k--','LineWidth',1);
grid on;
xlabel('Number of discretization points, N');
ylabel('Relative L2 error');
legend('40 fixed interior targets','Order 16: N^{-8}','Location','southwest');
title(sprintf('Torus Stokes %s: GL16 density, GL26 evaluation',layer));
% exportgraphics(gcf,fullfile(folder,sprintf('figure3_w16_generic_%s.png',layer)),'Resolution',200);

keyboard
