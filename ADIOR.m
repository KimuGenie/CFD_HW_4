%Alternating Direction Implicit Over-Relaxation method

function [u, residual] = ADIOR(u, dx, dy, imax, jmax, maxiter, tolerance, omega)

beta=dx/dy;
alpha = -2*(1+beta^2);
k=1;

Ai=zeros(imax-2);
for i=1:imax-3
    Ai(i, i) = alpha;
    Ai(i, i+1) = omega;
    Ai(i+1, i) = omega;
end
Ai(imax-2, imax-2) = alpha;

Ci=zeros(imax-2, 1);

Aj=zeros(jmax-2);
for j=1:jmax-3
    Aj(j, j) = alpha;
    Aj(j, j+1) = omega*beta^2;
    Aj(j+1, j) = omega*beta^2;
end
Aj(jmax-2, jmax-2) = alpha;

Cj=zeros(jmax-2, 1);

uprevi=u;

while k<=maxiter
    uprev=u;
    
    %x-sweep
    for j=2:jmax-1
        for i=2:imax-1
            Ci(i-1) = (1-omega)*alpha*u(i, j)-omega*beta^2*(uprev(i, j+1)+u(i, j-1));
        end
        Ci(1) = Ci(1)-omega*u(1, j);
        Ci(imax-2) = Ci(imax-2)-omega*u(imax, j);
        
        U=Ai\Ci; %Gauss elimination
        
        for i=2:imax-1
            uprevi(i, j) = U(i-1);
        end
    end
    
    for i=2:imax-1
        for j=2:jmax-1
            Cj(j-1) = (1-omega)*alpha*uprevi(i, j)-omega*(uprevi(i+1, j)+u(i-1, j));
        end
        Cj(1) = Cj(1)-omega*beta^2*u(i, 1);
        Cj(jmax-2) = Cj(imax-2)-omega*beta^2*u(i, jmax);
        
        U=Aj\Cj; %Gauss elimination
        
        for j=2:jmax-1
            u(i, j) = U(j-1);
        end
    end
    
    residual(k)=abs(mean(u-uprev, 'all'));
    if residual(k)<=tolerance
        break
    end
    k=k+1;
end