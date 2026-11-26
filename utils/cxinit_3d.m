function [center,cxloc,cxploc,cxpploc,fit]=cxinit_3d(T,x,Xu,Xv,Xuu,Xuv,Xvv)
% One GL16 patch: nodal position and local first/second derivatives (3 x 256 each, u fastest)
%   -> center and Legendre coefficients (P_0..P_15 x P_0..P_15) of x-center, and of the fitted derivatives.
% cxloc(16,16,3), cxploc(16,16,3,2), cxpploc(16,16,3,3); fit from compute_geometry_fit.
density_map=compute_legendre_values(15,T).'\eye(16);
center=mean(x,2);
fields={x-center,Xu,Xv,Xuu,Xuv,Xvv};
coeff=zeros(16,16,3,6);
for f=1:6
  for coordinate=1:3
    coeff(:,:,coordinate,f)=density_map.'*reshape(fields{f}(coordinate,:),16,16)*density_map;
  end
end
cxloc=coeff(:,:,:,1);
[cxploc,cxpploc,fit]=compute_geometry_fit(coeff(:,:,:,2:3),coeff(:,:,:,4:6),1.5,struct());
end
