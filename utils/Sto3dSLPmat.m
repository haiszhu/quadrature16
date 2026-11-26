function [A, T, Tx, Ty, Tz] = Sto3dSLPmat(t,s,if_self)
% dense Stokes SLP Nystrom eval matrix from sources to targets
%
% A = Sto3dSLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes
% Outputs:
%  A - (3*M*3*N) matrix getting potentials from density values
%
% NB: fills 3*M*3*N matrices, so assumes small cases. 
%

if nargin < 3, if_self = 0; end

x=t.x; y=s.x; M = size(x,2);
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
d1d1 = d1.^2; d2d2 = d2.^2; d3d3 = d3.^2;
r2 = d1d1+d2d2+d3d3;
ir = 1./sqrt(r2);   % 1/r matrix, size M*N
if if_self, ir(diagind(ir)) = 0; end
ir3w = ir.^3.*s.w(:)';
irw = ir.*s.w(:)';

d1d2 = d1.*d2; d1d3 = d1.*d3; d2d3 = d2.*d3;
d1d1ir3w = d1d1.*ir3w; d1d2ir3w = d1d2.*ir3w; d1d3ir3w = d1d3.*ir3w;
d2d2ir3w = d2d2.*ir3w; d2d3ir3w = d2d3.*ir3w; d3d3ir3w = d3d3.*ir3w;
Across = [d1d1ir3w,    d1d2ir3w,    d1d3ir3w;...
          d1d2ir3w,    d2d2ir3w,    d2d3ir3w;...
          d1d3ir3w,    d2d3ir3w,    d3d3ir3w];
A = 1/(8*pi)*(kron(eye(3),irw) + Across);   

if nargout == 2 % want traction, need t.nx
  d1nx1 = d1.*t.nx(1,:)'; % (d.nx)/r^2\
  d2nx2 = d2.*t.nx(2,:)';
  d3nx3 = d3.*t.nx(3,:)';
  dnxir2 = (d1nx1+d2nx2+d3nx3)./r2;
  if if_self, dnxir2(diagind(dnxir2)) = 0; end
  T = -3/(4*pi)*kron(ones(3),dnxir2).*Across;
end
if nargout == 5 % want traction x, y, z parts, multiplied by t.nx outside
  ir2 = ir.^2;
  d1ir2 = d1.*ir2; d2ir2 = d2.*ir2; d3ir2 = d3.*ir2;
  if if_self, d1ir2(diagind(d1ir2)) = 0; d2ir2(diagind(d2ir2)) = 0; d3ir2(diagind(d3ir2)) = 0; end
  Tx = -3/(4*pi)*kron(ones(3),d1ir2).*Across;
  Ty = -3/(4*pi)*kron(ones(3),d2ir2).*Across;
  Tz = -3/(4*pi)*kron(ones(3),d3ir2).*Across;
  T = kron([1;1;1],t.nx(1,:)').*Tx + kron([1;1;1],t.nx(2,:)').*Ty + kron([1;1;1],t.nx(3,:)').*Tz;
end

end

function i = diagind(A)
% function i = diagind(A)
%
% return diagonal indices of a square matrix, useful for changing a diagonal
% in O(N) effort, rather than O(N^2) if add a matrix to A using matlab diag()
%
% barnett 2/6/08

N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);
end