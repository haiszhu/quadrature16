function [x,nx,sp,dp,dt,rts,rps,rtts,rtps,rpps] = chestnutparam2(a,th,u,v,o,x,nx,sp,dp,dt,rts,rps,rtts,rtps,rpps)
% 2nd version, so old code does not break... need to fix creat panel
% routine 
%
% a = 1/2; 
% chestnutparam(u,v)
% rx = (...)*cos(u);
% ry = (...)*sin(u);
% rz = ...*cos(v);

% % th \in [pi/3, 2*pi/3];
% th = @(u) pi/6*sin(u) + pi/2; % need to be careful here

r = @(u,v)  [(a+sin(v/2).*cos(th*(v/(2*pi)-1/2))).*cos(u);...
             (a+sin(v/2).*cos(th*(v/(2*pi)-1/2))).*sin(u);...
             sin(v/2).*sin(th*(v/(2*pi)-1/2))]; % t,p rows
rt = @(u,v) [(1/2*cos(v/2).*cos(th*(v/(2*pi)-1/2)) - th/(2*pi)*sin(v/2).*sin(th*(v/(2*pi)-1/2))).*cos(u); ...
             (1/2*cos(v/2).*cos(th*(v/(2*pi)-1/2)) - th/(2*pi)*sin(v/2).*sin(th*(v/(2*pi)-1/2))).*sin(u);...
             1/2*cos(v/2).*sin(th*(v/(2*pi)-1/2)) + th/(2*pi)*sin(v/2).*cos(th*(v/(2*pi)-1/2))]; 
rp = @(u,v) [-(a+sin(v/2).*cos(th*(v/(2*pi)-1/2))).*sin(u);...
             (a+sin(v/2).*cos(th*(v/(2*pi)-1/2))).*cos(u);...
             0*u];

if isfield(o,'recenter') && o.recenter
  x = [sin(v/2).*cos(th*(v/(2*pi)-1/2)).*cos(u);
      sin(v/2).*cos(th*(v/(2*pi)-1/2)).*sin(u);
      sin(v/2).*sin(th*(v/(2*pi)-1/2))];
  x(1,:) = x(1,:)-2*a*sin(u/2).^2;
  x(2,:) = x(2,:)+a*sin(u);
else
  x = r(u,v);
end
rts = rt(u,v); rps = rp(u,v);
dt = sqrt(sum(rts.^2,1));
dp = sqrt(sum(rps.^2,1));    % could reuse this data
nx =  cross(rps, rts);       % outward normal
sp = sqrt(sum(nx.^2,1)); % speeds
nx = nx.*repmat(1./sp, [3 1]);
if nargout>7
  % Analytic Hessians avoid differentiating nearly constant thin-patch tangents.
  rate = th/(2*pi); phase = rate*v-th/2;
  rho_v = .5*cos(v/2).*cos(phase)-rate*sin(v/2).*sin(phase);
  rho_vv = -(.25+rate^2)*sin(v/2).*cos(phase)-rate*cos(v/2).*sin(phase);
  z_vv = -(.25+rate^2)*sin(v/2).*sin(phase)+rate*cos(v/2).*cos(phase);
  rtts = [rho_vv.*cos(u);rho_vv.*sin(u);z_vv];
  rtps = [-rho_v.*sin(u);rho_v.*cos(u);0*v];
  rpps = [-rps(2,:);rps(1,:);0*v];
end
end
