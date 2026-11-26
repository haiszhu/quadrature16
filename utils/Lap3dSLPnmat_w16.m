function [Anx,Any,Anz] = Lap3dSLPnmat_w16(t,s)
% LAP3DSLPNMAT_W16.  dense Laplace SLP target-gradient Nystrom eval matrices from sources to targets
%
% [Anx,Any,Anz] = Lap3dSLPnmat_w16(t,s) variant that requires parametric space info
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
%  Anx, Any, Anz - (M*N) matrices getting the x, y, z target derivatives from density values
%      (Lap3dSLPmat's Anx, Any, Anz); the caller contracts them with the target normals.
%
% NB: fills M*N matrices, so assumes small cases. 

%
M = size(t.x,2); N = size(s.x,2);
if M == 0
  Anx = zeros(0,N); Any = Anx; Anz = Anx;
  return
end
assert(N == 256, 'Lap3dSLPnmat_w16 requires one GL16 source panel.');
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
WSn = zeros(M,N,3);
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
tnx = zeros(3,M);
work = [];
ier = zeros(M,1);
[~,~,WSn,~,t0,h,work,ier] = quadSD_w16_mex(M,t.x,tnx,N,s.x,s.nx,[],[],WSn,[],t0,h,s.w,s.center,s.C_X,s.C_Xi,s.C_Xij,s.C_X_fit,cols,target_idx,t.otherside_flag,t.onsurface_comp,ib,xb,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,0,0,1,0,work,ier);

%
if any(ier ~= 0), error('Lap3dSLPnmat_w16:quadrature','Quadrature failed for one or more targets.'); end
Anx = 0.5*WSn(:,:,1); Any = 0.5*WSn(:,:,2); Anz = 0.5*WSn(:,:,3);

end