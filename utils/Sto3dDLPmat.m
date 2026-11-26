function [A,An] = Sto3dDLPmat(t,s,if_self)
% dense Stokes DLP Nystrom eval matrix from sources to targets
%
% A = Sto3dDLPmat(t,s)
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
% modified from SLP traction, add DLP traction later...
% DLP traction, modifed from Bobbie's Sto3dDLPmat, https://github.com/bobbielf2/ZetaTrap3D_Unified

if(nargin==0), testSto3dDLPn_mex(); return; end % this is to check Stokes DLPn fortran mex matches matlab output

if nargin < 3, if_self = 0; end

x=t.x; y=s.x; M = size(x,2);
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
d1d1 = d1.^2; d2d2 = d2.^2; d3d3 = d3.^2;
r2 = d1d1+d2d2+d3d3;
ir = 1./sqrt(r2);   % 1/r matrix, size M*N
if if_self, ir(diagind(ir)) = 0; end
ir2 = ir.^2;
ir3w = ir.*ir2.*s.w(:)';

d1d2 = d1.*d2; d1d3 = d1.*d3; d2d3 = d2.*d3;
d1d1ir3w = d1d1.*ir3w; d1d2ir3w = d1d2.*ir3w; d1d3ir3w = d1d3.*ir3w;
d2d2ir3w = d2d2.*ir3w; d2d3ir3w = d2d3.*ir3w; d3d3ir3w = d3d3.*ir3w;
Across = [d1d1ir3w,    d1d2ir3w,    d1d3ir3w;...
          d1d2ir3w,    d2d2ir3w,    d2d3ir3w;...
          d1d3ir3w,    d2d3ir3w,    d3d3ir3w];

dnxir2 = (d1.*s.nx(1,:)+d2.*s.nx(2,:)+d3.*s.nx(3,:)).*ir2;
if if_self, dnxir2(diagind(dnxir2)) = 0; end
A = 3/(4*pi)*kron(ones(3),dnxir2).*Across;

if nargout > 1 % 
  mu = 1;
  % tensor product d x d
  dxd = [d1d1,  d1d2,  d1d3; 
         d1d2,  d2d2,  d2d3;
         d1d3,  d2d3,  d3d3];
  ir5 = ir3w.*ir2;            % 1/r^5 * w, weights included
  dnyir5 = (d1.*s.nx(1,:)+d2.*s.nx(2,:)+d3.*s.nx(3,:)).*ir5;  % (d.ny)/r^5 * w
  nyir3 = {s.nx(1,:).*ir3w;s.nx(2,:).*ir3w;s.nx(3,:).*ir3w};
  dnx = d1.*t.nx(1,:)'+d2.*t.nx(2,:)'+d3.*t.nx(3,:)'; % d.nx
  dnxir5 = dnx.*ir5;
  ddnxir5 = {d1.*dnxir5;d2.*dnxir5;d3.*dnxir5}; % (d.nx)/r^5 * w
  ddnyir5 = {d1.*dnyir5;d2.*dnyir5;d3.*dnyir5}; % (d.ny)/r^5 * w
  nxnyir5 = (t.nx(1,:)'.*s.nx(1,:)+t.nx(2,:)'.*s.nx(2,:)+t.nx(3,:)'.*s.nx(3,:)).*ir5; % (nx.ny)/r^5 * w
  dnxdnyir5 = dnx.*dnyir5; % (d.nx)(d.ny)/r^5 * w
  An = kron((3*mu/4/pi)*ones(3),nxnyir5-10*dnxdnyir5.*ir2).*dxd; % (3(nx.ny)/r^5-30(d.nx)(d.ny)/r^7)*(d x d), incl prefac
  An = An + kron((3*mu/4/pi)*eye(3),dnxdnyir5); % + 3(d.nx)(d.ny)/r^5 * w * I
  % tensor prod terms
  TP = cell(3);
  for i = 1:3
    for j = 1:3
      TP{i,j} = 2*t.nx(i,:)'.*nyir3{j} + ...	 % 2*(nx)x(ny)/r^3 * w
                3*s.nx(i,:).*ddnxir5{j} + ...  % 3*(d.nx)*[(ny)x(d)]/r^5 * w
                3*ddnyir5{i}.*t.nx(j,:)';      % 3*(d.ny)*[(d)x(nx)]/r^5 * w
    end
  end
  An = An + (mu/4/pi)*cell2mat(TP); % DLP traction, incl prefac
  %keyboard
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

function testSto3dDLPn_mex()

addpath('../')
clear
setup()
qntype = 'T'; % or 'V'; tensor product or vioreanu-rokhlin nodes
order = 12; % G-S node order, or vioreanu nodes order
type = 'torus', so.a = 1; so.b = 0.5; o = []; o.fquad = false; o.p = order;
so.mp = 2*6; so.np = 2*10;
[pan0 N] = create_panels(type,so,o);
side = 'i';

% triangular or quadrilateral panels
pan = [];
for k = 1:numel(pan0)
  % original coordinate: quadrature and density sampling
  [sl,su] = get_high_order_quad(pan0{k},order,qntype);
  pan{k} = sl;
  h_dim = order^2;
end

k = so.mp;
pank = pan{k}; 
direc = 4*pank.nx*sqrt(sum(pank.w(:)));  if side == 'i', direc = -direc; end
t.x = [pank.x + 1/20*direc, pank.x + 1/10*direc, pank.x + 1/5*direc];
t.nx = rand(size(t.x));

% matlab
nt = numel(t.x(1,:)); ns = numel(pank.x(1,:));
[~,An] = Sto3dDLPmat(t,pank);
[An11 An12 An13 An22 An23 An33 An21 An31 An32] = Sto3dDLPnmat_mex(t,pank);


diff = abs(An - [An11 An12 An13;An21 An22 An23;An31 An32 An33]);
figure(1),clf,imagesc(log10(diff))

keyboard

end

function [sl,su] = get_high_order_quad(pan,order,qntype)

  [x,w]=gauss(order); [x1 x2] = meshgrid(x); ww = w(:)*w(:)'; 
  uvs = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
  wts = ww(:)';
  
  % 
  [rs, nrs, sps] = pan.chart(uvs); % these should be interpolated, update 'panel_smooth_quadk.m'
  ws = sps.*wts*pan.spfac;
  sl.x = rs; sl.nx = nrs; sl.w = ws; sl.bdchart = pan.bdchart; 
  su = [];

end