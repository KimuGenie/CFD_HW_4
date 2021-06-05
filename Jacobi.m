%Jacobi iteration method

function [u, residual] = Jacobi(u, dx, dy, imax, jmax, maxiter, tolerance)

beta=dx/dy;
k=1;

while k<=maxiter
    uprev=u;
    for j=2:jmax-1
        for i=2:imax-1
            u(i, j) = (uprev(i+1, j)+uprev(i-1, j)+beta^2*(uprev(i, j+1)+uprev(i, j-1)))/(2*(1+beta^2));
        end
    end
    
    residual(k)=mean(abs(u-uprev), 'all');
    if residual(k)<=tolerance
        break
    end
    k=k+1;
end

end