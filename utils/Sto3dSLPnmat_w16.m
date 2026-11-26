function [Tx,Ty,Tz] = Sto3dSLPnmat_w16(t,s)
% STO3DSLPNMAT_W16.  dense Stokes SLP traction Nystrom eval matrices from sources to targets
%
% [Tx,Ty,Tz] = Sto3dSLPnmat_w16(t,s) variant that requires parametric space info
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
%      x (3*N) nodes
%      onsurface_comp (1) scalar logical; true for boundary evaluation
%      otherside_flag (1*M) logical; true for targets on the opposite side relative to s
%
% Outputs:
%  Tx, Ty, Tz - (3M*3N) matrices, blocks [u1;u2;u3] x [f1;f2;f3] as Sto3dSLPmat, getting the x, y, z parts of the SLP
%      traction -3/(4*pi)*((x-y).n)(x-y)(x-y)'/|x-y|^5 from density values (Sto3dSLPmat's Tx, Ty, Tz, viscosity 1);
%      the caller contracts them with the target normals, T = kron([1;1;1],t.nx(1,:)').*Tx + ... .
%      Tk is quadSD_stok_w16_mex's SLPn with the target normal tnx = e_k: exact parts off surface and in the mid; on surface
%      (h=0 rows, the direct driver at the chart normal n_r0) the normal part (e_k.n_r0)*T(n_r0), so the contraction with
%      the target normal t.nx = +-n_r0 is exact there too.
%
% NB: fills 3M*3N matrices, so assumes small cases.

%
M = size(t.x,2); N = size(s.x,2);
if M == 0
  Tx = zeros(0,3*N); Ty = Tx; Tz = Tx;
  return
end
assert(N == 256, 'Sto3dSLPnmat_w16 requires one GL16 source panel.');
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
T = cell(1,3);
for k = 1:3
  tnx = zeros(3,M); tnx(k,:) = 1;   % target normal e_k
  WSn = zeros(M,N,3,3);
  ier = zeros(M,1);
  [~,~,WSn,~,~,~,~,ier] = quadSD_stok_w16_mex(M,t.x,tnx,N,s.x,s.nx,[],[],WSn,[],t0,h,s.w,s.center,s.C_X,s.C_Xi,s.C_Xij,s.C_X_fit,cols,target_idx,t.otherside_flag,t.onsurface_comp,ib,xb,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,0,0,1,0,work,ier);
  if any(ier ~= 0), error('Sto3dSLPnmat_w16:quadrature','Quadrature failed for one or more targets.'); end
  T{k} = reshape(permute(WSn,[1 3 2 4]),3*M,3*N);   % WSn(i,j,a,b) -> T{k}(i+(a-1)*M,j+(b-1)*N)
end
[Tx,Ty,Tz] = T{:};

end
