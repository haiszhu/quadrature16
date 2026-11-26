function delta = chestnutdisplacement(a,th,v,u0,v0,du,orientation,orientation0,delta)
% Source minus target, retaining the supplied azimuth difference du before position rounding.
% v: source column; u0,v0: target rows; du: source-by-target; delta: 3-by-source-by-target.
% orientation0 is a scalar or a row containing one reflection sign per target.
rho = sin(v/2).*cos(th*(v/(2*pi)-.5));
z = sin(v/2).*sin(th*(v/(2*pi)-.5));
rho0 = sin(v0/2).*cos(th*(v0/(2*pi)-.5));
z0 = sin(v0/2).*sin(th*(v0/(2*pi)-.5));
radial = (rho-rho0)-2*(a+rho).*sin(du/2).^2;
azimuthal = (a+rho).*sin(du);
delta(1,:,:) = reshape(cos(u0).*radial-sin(u0).*azimuthal,1,numel(v),[]);
delta(2,:,:) = reshape(sin(u0).*radial+cos(u0).*azimuthal,1,numel(v),[]);
delta(3,:,:) = reshape(orientation*z-orientation0.*z0,1,numel(v),[]);
end
