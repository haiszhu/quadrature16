function [r, nr, w] = panel_smooth_quadk(chart,psiz,p)   % setup G-L tensor native quadr
[x, w, D] = compute_gauss_legendre(p);   % Gauss-Legendre nodes, weights, differentiation matrix on [-1,1]
[x1 x2] = meshgrid(x); xx = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
ww = w(:)*w(:)'; ww = ww(:)';                     % 1*p^2 weights

  % only use pan{m}.x
  r = chart(xx);  % get 3
  % compute pan{m}.nx & sp from r
  tp_ind = reshape(1:p^2,[p p]);
  xts = nan(1,p^2); yts = xts; zts = xts; % partial t (small circle)
  for k=1:p
    xk = r(1,tp_ind(:,k)); yk = r(2,tp_ind(:,k)); zk = r(3,tp_ind(:,k));
    xts(tp_ind(:,k)) = 1/psiz(2)*D*xk(:); yts(tp_ind(:,k)) = 1/psiz(2)*D*yk(:); zts(tp_ind(:,k)) = 1/psiz(2)*D*zk(:);
  end
  rts = [xts(:),yts(:),zts(:)]';
  xps = nan(1,p^2); yps = xps; zps = xps; % partial p (big circle)
  for k=1:p
    x = r(1,tp_ind(k,:)); y = r(2,tp_ind(k,:)); z = r(3,tp_ind(k,:));
    xps(tp_ind(k,:)) = 1/psiz(1)*D*x(:); yps(tp_ind(k,:)) = 1/psiz(1)*D*y(:); zps(tp_ind(k,:)) = 1/psiz(1)*D*z(:);
  end
  rps = [xps(:),yps(:),zps(:)]';
  nr =  cross(rps, rts);       % outward normal
  sp = sqrt(sum(nr.^2,1)); % speeds
  nr = nr.*repmat(1./sp, [3 1]);
  spfac = psiz(1)*psiz(2);
  w = sp.*ww*spfac;
  t = xx;          % same for all panels, for now
  
