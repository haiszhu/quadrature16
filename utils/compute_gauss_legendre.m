function [x,w,D] = compute_gauss_legendre(n)
% COMPUTE_GAUSS_LEGENDRE  n-point Gauss-Legendre rule on [-1,1], n >= 2: nodes x_k (ascending) and weights w_k, columns.
% DLMF 3.5(v): the x_k are the zeros of p_n = P_n, and the rule is exact for polynomials of degree <= 2n-1.
%   Symmetry: x_{n+1-k} = -x_k, w_{n+1-k} = w_k; for odd n the middle node is 0.
%   Nodes k <= n/2: Newton's method, x <- x-p_n(x)/p_n'(x), from x_k = -cos((k-1/4)*pi/(n+1/2)), until no node
%   moves by more than eps; p_{n-1}, p_n by DLMF 18.9.2 as in compute_legendre_values, and
%   p_n'(x) = n*(p_{n-1}(x)-x*p_n(x))/(1-x^2) (DLMF 14.10.5, mu = 0).
%   Weights: by DLMF (3.5.20_2), sum_k p_l(x_k)*p_m(x_k)*w_k = h_m*delta_lm for l, m <= n-1, with h_m = 2/(2m+1)
%   (DLMF Table 18.3.1), so 1/w_k = sum_{m=0}^{n-1} (m+1/2)*p_m(x_k)^2.
% D (when requested): D(i,j) = derivative at x_i of the polynomial interpolating at the nodes, D*f = f' at the nodes.
%   With P(k,m+1) = p_m(x_k), (3.5.20_2) reads P.'*diag(w)*P = diag(h), so the Legendre coefficients of f are
%   diag(1/h)*P.'*diag(w)*f and D = P'*diag(1/h)*P.'*diag(w), P'(k,m+1) = p_m'(x_k).
a_n=((0:n)+1)./(2*(0:n)+1);   % index m+1 holds a_m, b_m, c_m
b_n=zeros(1,n+1);
c_n=(0:n)./(2*(0:n)+1);
a_n_inv=1./a_n;
n_half=floor(n/2);
n_w=ceil(n/2);                % nodes x <= 0
x=zeros(n,1);
x(1:n_half)=-cos(((1:n_half).'-0.25)*pi/(n+0.5));
for iteration=1:10
  x_half=x(1:n_half);
  p_mm1=zeros(n_half,1);      % p_{-1}
  p_m=ones(n_half,1);         % p_0
  for m=0:n-1
    p_mp1=((x_half-b_n(m+1)).*p_m-c_n(m+1)*p_mm1)*a_n_inv(m+1);
    p_mm1=p_m;
    p_m=p_mp1;
  end
  % p_m = p_n, p_mm1 = p_{n-1}
  x_new=x_half-p_m./(n*(p_mm1-x_half.*p_m)./((1-x_half).*(1+x_half)));
  step=max(abs(x_new-x_half));
  x(1:n_half)=x_new;
  if step<=eps, break; end
end
if n_w>n_half, x(n_w)=0; end
p_mm1=zeros(n_w,1);
p_m=ones(n_w,1);
w=zeros(n,1);
for m=0:n-1
  w(1:n_w)=w(1:n_w)+(m+0.5)*p_m.^2;
  p_mp1=((x(1:n_w)-b_n(m+1)).*p_m-c_n(m+1)*p_mm1)*a_n_inv(m+1);
  p_mm1=p_m;
  p_m=p_mp1;
end
w(1:n_w)=1./w(1:n_w);
x(n+1-(1:n_half))=-x(1:n_half);
w(n+1-(1:n_half))=w(1:n_half);
if nargout<3, return; end
[p,dp]=compute_legendre_values(n-1,x);
D=(dp.*((0:n-1)+0.5))*(p.'.*w.');
end
