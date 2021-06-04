%Line Successive Over-Relaxation method x-sweep

function [u, residual] = LSORx(u, dx, dy, imax, jmax, maxiter, tolerance, omega)

beta=dx/dy;
alpha = -2*(1+beta^2);
k=1;

A=zeros(imax-2);
for i=1:imax-3
    A(i, i) = alpha;
    A(i, i+1) = omega;
    A(i+1, i) = omega;
end
A(imax-2, imax-2) = alpha;

C=zeros(imax-2, 1);

while k<=maxiter
    uprev=u;
    for j=2:jmax-1
        for i=2:imax-1
            C(i-1) = (1-omega)*alpha*u(i, j)-omega*beta^2*(uprev(i, j+1)+u(i, j-1));
        end
        C(1) = C(1)-omega*u(1, j);
        C(imax-2) = C(imax-2)-omega*u(imax, j);
        
        U=A\C; %Gauss elimination
        
        for i=2:imax-1
            u(i, j) = U(i-1);
        end
    end
    residual(k)=mean(abs(u-uprev), 'all');
    if residual(k)<=tolerance
        break
    end
    k=k+1;
end

end