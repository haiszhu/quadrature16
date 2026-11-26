function An = Sto3dDLPnmat_w16(t,s)
% STO3DDLPNMAT_W16.  dense Stokes DLP traction Nystrom eval matrix from sources to targets
%
% An = Sto3dDLPnmat_w16(t,s) variant that requires parametric space info
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
%      center (3*1) center of geometric expansion
%      C_X (16*16*3) Legendre coefficients of x-center
%      C_Xi (16*16*3*2) first derivative coefficients: u, v
%      C_Xij (16*16*3*3) second derivative coefficients: uu, uv, vv
%      C_X_fit (16*16*3) derivative-fitted position coefficients for projection
%
%  t - target struct with fields:
%      x (3*M) nodes, nx (3*M) unit target normals
%      onsurface_comp (1) scalar logical; true for boundary evaluation
%      otherside_flag (1*M) logical; true for targets on the opposite side relative to s
%
% Outputs:
%  An - (3M*3N) matrix getting the traction of the DLP at the target normal t.nx from density values, blocks
%      [t1;t2;t3] x [f1;f2;f3] as Sto3dDLPmat's second output (viscosity 1, source normal n = s.nx). It is -WDn of
%      quadSD_stok_w16_mex: off surface and in the mid the kernel at t.nx; on surface (h=0 rows) the finite part at
%      the chart normal n_r0 times t.nx.n_r0, exact for t.nx = +-n_r0 (no jump term).
%
% NB: fills 3M*3N matrices, so assumes small cases.

%
M = size(t.x,2); N = size(s.x,2);
if M == 0
  An = zeros(0,3*N);
  return
end
assert(N == 256, 'Sto3dDLPnmat_w16 requires one GL16 source panel.');
assert(isscalar(t.onsurface_comp));
if ~isfield(t,'otherside_flag')
  t.otherside_flag = false(1,M);
end
assert(numel(t.otherside_flag) == M);
t0 = nan(2,M);
h  = nan(M,1);
if isfield(t,'t0'), t0 = t.t0;  end
if isfield(t,'h'),  h = t.h(:); end
if_h0 = t.onsurface_comp & ~t.otherside_flag;
WDn = zeros(M,N,3,3);
cols = 1:N;
target_idx = zeros(1,M);
if t.onsurface_comp
  [~,idx] = ismember(t.x.',s.x.','rows');
  target_idx = idx.';
  target_idx(logical(t.otherside_flag)) = 0;
end

%
ib = 0; xb = []; cols_typeb = [];
cxb = []; cxpb = []; cxppb = [];
tcen = []; tsca = []; tau0 = nan(2,M);

%
work = [];
ier = zeros(M,1);
[~,~,~,WDn,t0,h,work,ier] = quadSD_stok_w16_mex(M,t.x,t.nx,N,s.x,s.nx,[],[],[],WDn,t0,h,s.w,s.center,s.C_X,s.C_Xi,s.C_Xij,s.C_X_fit,cols,target_idx,t.otherside_flag,t.onsurface_comp,ib,xb,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,0,0,0,1,work,ier);

%
if any(ier ~= 0), error('Sto3dDLPnmat_w16:quadrature','Quadrature failed for one or more targets.'); end
An = -reshape(permute(WDn,[1 3 2 4]),3*M,3*N);   % WDn(i,j,ic,jc) -> -An(i+(ic-1)*M,j+(jc-1)*N)

end
