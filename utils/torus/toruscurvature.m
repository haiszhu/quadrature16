function [hx] = toruscurvature(a,b,p,t,o)
% TORUSPARAM - return coords, unit normal and speed at params (phi,theta)
%
% [x nx sp dp dt] = torusparam(a,b,p,t).
%  Inputs:
%  p,t are phi and theta, and should be row vectors if multiple points wanted
%  a is major radius
%  b is minor radius.
%    If b is a cell array of 3 function handles, is interpreted as b(t,p),
%    a minor (poloidal) radius function over theta and phi each in [0,2pi),
%    and b_t(t,p), b_p(t,p) its correct partials. Or, see below option.
%
% Geometry:
%  1st coord (phi) goes CCW (viewed from above) around big circle, toroidal
%  2nd coord (theta) takes up around little circle, poloidal
% (dp,dt,outwardsnormal) form a RH coord sys.
%
% sp = dp * dt
% dt is magnitude of dr/dt vector; dp mag of dr/dp vector.
%
% [x nx sp dp dt] = torusparam(a,b,p,t,o) allows opts eg
% o.f, o.fp, o.ft handles to modulation func (poloidal radius vs theta,phi),
% and its partials, which must be correct.

% adapted from Barnett 2013 except made outward normal dphi cross dtheta.
% expanded modulation interface 10/10/18.

if nargin>4 || iscell(b)       % wobbly torus
    if iscell(b)    % allow a new interface to wobbly modulation
    f = b{1}; ft = b{2}; fp = b{3}; ftt = b{4}; ftp = b{5}; fpp = b{6};
    else
    f = o.f; ft = o.ft; fp = o.fp;
    end
%     r = @(t,p) [(a+f(t,p).*cos(t)).*cos(p); (a+f(t,p).*cos(t)).*sin(p); f(t,p).*sin(t)]; % t,p rows
%     rt = @(t,p) [(-f(t,p).*sin(t) + ft(t,p).*cos(t)).*cos(p); (-f(t,p).*sin(t)+ft(t,p).*cos(t)).*sin(p); f(t,p).*cos(t)+ft(t,p).*sin(t)]; % partial
%     rp = @(t,p) [-(a+f(t,p).*cos(t)).*sin(p)+fp(t,p).*cos(t).*cos(p); (a+f(t,p).*cos(t)).*cos(p) + fp(t,p).*cos(t).*sin(p); fp(t,p).*sin(t)];
%     rtt = @(t,p) [(ftt(t,p).*cos(t)-2*ft(t,p).*sin(t)-f(t,p).*cos(t)).*cos(p); (ftt(t,p).*cos(t)-2*ft(t,p).*sin(t)-f(t,p).*cos(t)).*sin(p); ftt(t,p).*sin(t)+2*ft(t,p).*cos(t)-f(t,p).*sin(t)];
%     rtp = @(t,p) [ftp(t,p).*cos(t).*cos(p)-fp(t,p).*sin(t).*cos(p)-ft(t,p).*cos(t).*sin(p)+f(t,p).*sin(t).*sin(p);...
%                 ftp(t,p).*cos(t).*sin(p)-fp(t,p).*sin(t).*sin(p)+ft(t,p).*cos(t).*cos(p)-f(t,p).*sin(t).*cos(p);...
%                 ftp(t,p).*sin(t)+fp(t,p).*cos(t)];
%     rpp = @(t,p) [fpp(t,p).*cos(t).*cos(p)-2*fp(t,p).*cos(t).*sin(p)-(a+f(t,p).*cos(t)).*cos(p);...
%                 fpp(t,p).*cos(t).*sin(p)+2*fp(t,p).*cos(t).*cos(p)-(a+f(t,p).*cos(t)).*sin(p);...
%                 fpp(t,p).*sin(t)];
%     fnormal = @(t,p) [f(t,p).*fp(t,p).*sin(p)+(a+f(t,p).*cos(t)).*(ft(t,p).*sin(t)+f(t,p).*cos(t)).*cos(p);...
%                    -f(t,p).*fp(t,p).*cos(p)+(a+f(t,p).*cos(t)).*(ft(t,p).*sin(t)+f(t,p).*cos(t)).*sin(p);...
%                    -(a+f(t,p).*cos(t)).*(ft(t,p).*cos(t)-f(t,p).*sin(t))]./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));

    f2N = @(t,p) (-2*f(t,p).*fp(t,p).^2.*cos(t)-(a+f(t,p).*cos(t)).^2.*(ft(t,p).*sin(t)+f(t,p).*cos(t))+fpp(t,p).*f(t,p).*(a+f(t,p).*cos(t)))...
                 ./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));
    f2M = @(t,p) (f(t,p).^2.*fp(t,p).*sin(t)-f(t,p).*ft(t,p).*fp(t,p).*cos(t)+(ftp(t,p).*f(t,p)-ft(t,p).*fp(t,p)).*(a+f(t,p).*cos(t)))...
                 ./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));
    f2L = @(t,p) (a+f(t,p).*cos(t)).*(f(t,p).*ftt(t,p)-2*ft(t,p).^2-f(t,p).^2)./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));

    f1G = @(t,p) fp(t,p).^2+(a+f(t,p).*cos(t)).^2;
    f1F = @(t,p) ft(t,p).*fp(t,p);
    f1E = @(t,p) ft(t,p).^2+f(t,p).^2;

    fH = @(t,p) (f1E(t,p).*f2N(t,p)-2*f1F(t,p).*f2M(t,p)+f1G(t,p).*f2L(t,p))./(2*(f1E(t,p).*f1G(t,p)-f1F(t,p).^2));
else % param...
%     f = @(t,p) [(a+b*cos(t)).*cos(p); (a+b*cos(t)).*sin(p); b*sin(t)]; % t,p rows
%     ft = @(t,p) [-b*sin(t).*cos(p); -b*sin(t).*sin(p); b*cos(t)]; % partial
%     fp = @(t,p) [-(a+b*cos(t)).*sin(p); (a+b*cos(t)).*cos(p); 0*p];
%     ftt = @(t,p) [-b*cos(t).*cos(p); -b*cos(t).*sin(p); -b*sin(t)];
%     fpp = @(t,p) [-(a+b*cos(t)).*cos(p); -(a+b*cos(t)).*sin(p); 0*p];
%     ftp = @(t,p) [b*sin(t).*sin(p); -b*sin(t).*cos(p); 0*p];
%     f2N = @(t,p) (-2*f(t,p).*fp(t,p).^2.*cos(t)-(a+f(t,p).*cos(t)).^2.*(ft(t,p).*sin(t)+f(t,p).*cos(t))+fpp(t,p).*f(t,p).*(a+f(t,p).*cos(t)))...
%                  ./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));
%     f2M = @(t,p) (f(t,p).^2.*fp(t,p).*sin(t)-f(t,p).*ft(t,p).*fp(t,p).*cos(t)+(ftp(t,p).*f(t,p)-ft(t,p).*fp(t,p)).*(a+f(t,p).*cos(t)))...
%                  ./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));
%     f2L = @(t,p) (a+f(t,p).*cos(t)).*(f(t,p).*ftt(t,p)-2*ft(t,p).^2-f(t,p).^2)./sqrt(f(t,p).^2.*fp(t,p).^2+(a+f(t,p).*cos(t)).^2.*(ft(t,p).^2+f(t,p).^2));
% 
%     f1G = @(t,p) fp(t,p).^2+(a+f(t,p).*cos(t)).^2;
%     f1F = @(t,p) ft(t,p).*fp(t,p);
%     f1E = @(t,p) ft(t,p).^2+f(t,p).^2;
% 
%     fH = @(t,p) (f1E(t,p).*f2N(t,p)-2*f1F(t,p).*f2M(t,p)+f1G(t,p).*f2L(t,p))./(2*(f1E(t,p).*f1G(t,p)-f1F(t,p).^2));
    f2N = @(t,p) -b*(a+b*cos(t)).^2.*cos(t)./abs((a+b*cos(t))*b);
    f2M = @(t,p) 0*t;
    f2L = @(t,p) -b^2*(a+b*cos(t))./abs((a+b*cos(t))*b);
  
    f1G = @(t,p) (a+b*cos(t)).^2;
    f1F = @(t,p) 0*t;
    f1E = @(t,p) b^2*ones(size(t));

    fH = @(t,p) (f1E(t,p).*f2N(t,p)-2*f1F(t,p).*f2M(t,p)+f1G(t,p).*f2L(t,p))./(2*(f1E(t,p).*f1G(t,p)-f1F(t,p).^2));
%     keyboard
end

hx = fH(t,p);
% rts = rt(t,p); rps = rp(t,p);
% dt = sqrt(sum(rts.^2,1));
% dp = sqrt(sum(rps.^2,1));    % could reuse this data
% nx =  cross(rps, rts);       % outward normal
% sp = sqrt(sum(nx.^2,1)); % speeds
% nx = nx.*repmat(1./sp, [3 1]);