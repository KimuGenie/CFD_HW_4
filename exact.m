function u = exact(H, L, dx, dy, imax, jmax)

u=zeros(imax, jmax);
u(:, 1)=100;
T1=100;
for j=2:jmax
    y=(j-1)*dy;
    for i=1:imax-1
        residual=1;
        n=1;
        x=(i-1)*dx;
        while residual>eps
            uprev=u(i, j);
            u(i, j)=u(i, j)+T1*2*(1-(-1)^(n))/(n*pi)*sinh((n*pi*(H-y))/L)/sinh((n*pi*H)/L)*sin((n*pi*x)/L);
            residual=abs(u(i, j)-uprev);
            n=n+2;
        end
    end
end

end