function [C_Xi,C_Xij,fit] = compute_geometry_fit(C_Xi,C_Xij,residual_factor,fit)
%COMPUTE_GEOMETRY_FIT  Consistent source geometry: the five derivative fields of one polynomial X_h.
% Input: C_Xi(16,16,3,2) (X_1,X_2) and C_Xij(16,16,3,3) (X_11,X_12,X_22), the GL16 Legendre coefficients of
% the supplied nodal derivatives of one source patch; residual_factor>=1 (e.g. 1.5).
% Inout: C_Xi, C_Xij, overwritten by the coefficients of d_i X_h and d_ij X_h; fit.L, fit.residual
% (largest scaled nodal misfit of the selected fit) and fit.full_residual (the same for L=30).
% X_h(z)=sum over (m,n) in S_L of a_mn*P_m(z1)*P_n(z2), S_L={0<=m,n<=15, 1<=m+n<=L} (a_00 affects no
% derivative). fit.C_X(16,16,3) returns these position coefficients with a_00=0.
% a minimizes sum over f in {1,2,11,12,22} and the GL16 nodes g of
% |d_f X_h(g_i,g_j)-X_f(g_i,g_j)|^2/scale_f^2, scale_f=max|X_f| over the nodes and xyz (an identically zero
% field, e.g. the Hessian of a flat patch, takes the largest scale of the other fields).
% L is the smallest of 15,...,30 whose largest scaled misfit is at most residual_factor times that of L=30.
% Only the supplied derivatives of the source patch are used: no positions, targets or neighbours.
% The independently interpolated fields are not derivatives of one function (d_j X_i ~= X_ij); outside
% [-1,1]^2 that mismatch and the noise of high mixed modes grow quickly. The returned fields are exact
% derivatives of X_h, with the mixed modes m+n>L removed. Call with no arguments for the DLP/SLP checks.

%% 1. Nodal data and the derivative operators of the Legendre basis at the GL16 nodes.
[V,D1,D2]=compute_legendre_values(15,compute_gauss_legendre(16));    % V(i,m+1)=P_m(g_i), D1 = P_m', D2 = P_m'' there
Dm=zeros(16);                                            % d/dz P_n = sum_m Dm(m+1,n+1)*P_m
for degree=1:15
    for ell=degree-1:-2:0
        Dm(ell+1,degree+1)=2*ell+1;
    end
end
% Rows: X_1, X_2, X_11, X_12, X_22 at the 256 nodes; columns: vec(a), a(m+1,n+1).
M=[kron(V,D1);kron(D1,V);kron(V,D2);kron(D1,D1);kron(D2,V)];
fields={C_Xi(:,:,:,1),C_Xi(:,:,:,2),C_Xij(:,:,:,1),C_Xij(:,:,:,2),C_Xij(:,:,:,3)};
data=zeros(5*256,3); field_scale=zeros(5,1);
for f=1:5
    rows=(f-1)*256+(1:256);
    for coordinate=1:3
        data(rows,coordinate)=reshape(V*fields{f}(:,:,coordinate)*V.',[],1);
    end
    field_scale(f)=max(abs(data(rows,:)),[],'all');
end
field_scale(field_scale==0)=max(field_scale);
row_scale=repelem(1./field_scale,256);
M=row_scale.*M; data=row_scale.*data;
column_scale=2.^floor(log2(vecnorm(M))); column_scale(column_scale==0)=1;
M=M./column_scale; % Powers of two balance thin-patch columns without rounding their entries.

%% 2. Least squares on S_L; select the smallest L within residual_factor of the full set.
[m,n]=ndgrid(0:15,0:15);
fit.full_residual=misfit(M,data,m+n>=1);
for L=15:30
    selected=m+n>=1 & m+n<=L;
    fit.residual=misfit(M,data,selected);
    if fit.residual<=residual_factor*fit.full_residual, break, end
end
fit.L=L;
a=zeros(256,3); a(selected(:),:)=(M(:,selected(:))\data)./column_scale(selected(:)).';
fit.C_X=reshape(a,16,16,3);

%% 3. Coefficients of the derivatives of X_h: d_1 -> Dm*A, d_2 -> A*Dm'.
for coordinate=1:3
    A=reshape(a(:,coordinate),16,16);
    C_Xi(:,:,coordinate,1)=Dm*A; C_Xi(:,:,coordinate,2)=A*Dm.';
    C_Xij(:,:,coordinate,1)=Dm*(Dm*A); C_Xij(:,:,coordinate,2)=Dm*A*Dm.'; C_Xij(:,:,coordinate,3)=(A*Dm.')*Dm.';
end
end


function value=misfit(M,data,selected)
% Largest scaled nodal misfit of the least-squares fit on the selected modes.
a=M(:,selected(:))\data;
value=max(abs(M(:,selected(:))*a-data),[],'all');
end
