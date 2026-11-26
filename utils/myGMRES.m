function [x,it,predicted_relres,true_relres] = myGMRES(A,b,n,m,tol,shift)
%MYGMRES  Helsing GMRES for (shift*I+A)*x=b; shift defaults to 1.
% A is a matrix or a function handle x -> A*x.
% Adapted from Johan Helsing's demo1, myGMRES and rotmat:
% https://www.maths.lth.se/na/staff/helsing/Tutor/demo1.txt
% The tolerance applies to the predicted residual; the true residual is
% computed separately. M is the maximum number of iterations, no restart.
if nargin<6, shift=1; end

V=zeros(n,m+1);
H=zeros(m);
cs=zeros(m,1);
sn=zeros(m,1);
bnrm2=norm(b);
V(:,1)=b/bnrm2;
s=bnrm2*eye(m+1,1);
for it=1:m
    it1=it+1;
    if isnumeric(A), w=A*V(:,it); else, w=A(V(:,it)); end
    for k=1:it
        H(k,it)=V(:,k)'*w;
        w=w-H(k,it)*V(:,k);
    end
    H(it,it)=H(it,it)+shift;
    wnrm2=norm(w);
    V(:,it1)=w/wnrm2;
    for k=1:it-1
        temp=cs(k)*H(k,it)+sn(k)*H(k+1,it);
        H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
        H(k,it)=temp;
    end
    % Helsing's rotmat, inlined to keep the solver in one function.
    aa=H(it,it); bb=wnrm2;
    if bb==0
        cs(it)=1; sn(it)=0;
    elseif abs(bb)>abs(aa)
        temp=aa/bb;
        sn(it)=1/sqrt(1+temp^2);
        cs(it)=temp*sn(it);
    else
        temp=bb/aa;
        cs(it)=1/sqrt(1+temp^2);
        sn(it)=temp*cs(it);
    end
    H(it,it)=cs(it)*H(it,it)+sn(it)*wnrm2;
    s(it1)=-sn(it)*s(it);
    s(it)=cs(it)*s(it);
    predicted_relres=abs(s(it1))/bnrm2;
    if predicted_relres<=tol || it==m
        y=triu(H(1:it,1:it))\s(1:it);
        x=fliplr(V(:,1:it))*flipud(y);
        if isnumeric(A), Ax=A*x; else, Ax=A(x); end
        true_relres=norm(shift*x+Ax-b)/bnrm2;
        if nargout<3
            fprintf('  myGMRES: %d iterations, predicted residual %.3e, true residual %.3e\n', ...
                it,predicted_relres,true_relres);
        end
        break
    end
end
end