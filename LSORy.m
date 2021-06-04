%Line Successive Over-Relaxation method y-sweep

function [u, residual] = LSORy(u, dx, dy, imax, jmax, maxiter, tolerance, omega)

beta=dx/dy;
alpha = -2*(1+beta^2);
k=1;

A=zeros(imax-2);
for j=1:jmax-3
    A(j, j) = alpha;
    A(j, j+1) = omega;
    A(j+1, j) = omega;
end
A(jmax-2, jmax-2) = alpha;

C=zeros(jmax-2, 1);

while k<=maxiter
    uprev=u;
    for i=2:imax-1
        for j=2:jmax-1
            C(j-1) = (1-omega)*alpha*u(i, j)-omega*beta^2*(uprev(i+1, j)+u(i-1, j));
        end
        C(1) = C(1)-omega*u(i, 1);
        C(jmax-2) = C(jmax-2)-omega*u(i, jmax);
        
        U=A\C; %Gauss elimination
        
        for j=2:jmax-1
            u(i, j) = U(j-1);
        end
    end
    residual(k)=mean(abs(u-uprev), 'all');
    if residual(k)<=tolerance
        break
    end
    k=k+1;
end

end