function [p,dp,d2p] = compute_legendre_values(deg,x,d)
% COMPUTE_LEGENDRE_VALUES  Legendre polynomials p_n = P_n and their derivatives at points, n = 0..deg (deg >= 1).
%   p(:,n+1) = p_n(x) at the points x (any shape, taken as a column): numel(x) x (deg+1).
%   dp  = d.*p_n'(x)    = d/dlambda p_n(x+lambda*d) at lambda = 0;  d scalar or one per point, default 1.
%   d2p = d.^2.*p_n''(x). Only the requested outputs are formed. Valid at x = +-1 and outside [-1,1].
% DLMF 18.9.2, x*p_n = a_n*p_{n+1}+b_n*p_n+c_n*p_{n-1}, a_n = (n+1)/(2n+1), b_n = 0, c_n = n/(2n+1),
% solved for p_{n+1} (a_n_inv = 1/a_n):
%   p_{n+1}   = ((x-b_n)*p_n-c_n*p_{n-1})*a_n_inv,                     p_0 = 1, p_{-1} = 0,
%   dp_{n+1}  = ((x-b_n)*dp_n+d*p_n-c_n*dp_{n-1})*a_n_inv,             (its derivative along d)
%   d2p_{n+1} = ((x-b_n)*d2p_n+2*d*dp_n-c_n*d2p_{n-1})*a_n_inv.
a_n=((0:deg)+1)./(2*(0:deg)+1);   % index n+1 holds a_n, b_n, c_n
b_n=zeros(1,deg+1);
c_n=(0:deg)./(2*(0:deg)+1);
a_n_inv=1./a_n;
x=x(:);
p=zeros(numel(x),deg+1);
p(:,1)=1;
p_nm1=zeros(numel(x),1);
for n=0:deg-1
  p(:,n+2)=((x-b_n(n+1)).*p(:,n+1)-c_n(n+1)*p_nm1)*a_n_inv(n+1);
  p_nm1=p(:,n+1);
end
if nargout<2, return; end
if nargin<3, d=1; end
d=d(:);
dp=zeros(numel(x),deg+1);
dp_nm1=zeros(numel(x),1);
for n=0:deg-1
  dp(:,n+2)=((x-b_n(n+1)).*dp(:,n+1)+d.*p(:,n+1)-c_n(n+1)*dp_nm1)*a_n_inv(n+1);
  dp_nm1=dp(:,n+1);
end
if nargout<3, return; end
d2p=zeros(numel(x),deg+1);
d2p_nm1=zeros(numel(x),1);
for n=0:deg-1
  d2p(:,n+2)=((x-b_n(n+1)).*d2p(:,n+1)+2*d.*dp(:,n+1)-c_n(n+1)*d2p_nm1)*a_n_inv(n+1);
  d2p_nm1=d2p(:,n+1);
end
end
