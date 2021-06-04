%Jacobi iteration method

function [u, residual] = Jacobi(u, dx, dy, imax, jmax, maxiter, tolerance)

beta=dx/dy;
k=1;
ucurr=u;

while k<=maxiter
    uprev=u;
    for j=2:jmax-1
        for i=2:imax-1
            ucurr(i, j) = (u(i+1, j)+u(i-1, j)+beta^2*(u(i, j+1)+u(i, j-1)))/(2*(1+beta^2));
        end
    end
    u=ucurr;
    residual(k)=mean(abs(u-uprev), 'all');
    if residual(k)<=tolerance
        break
    end
    k=k+1;
end

end