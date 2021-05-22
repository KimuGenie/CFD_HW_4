function u = Five_point(u, dx, dy, imax, jmax)

beta = dx/dy;
alpha = -2*(1+beta^2);
size = (imax-2)*(jmax-2);
A = zeros(size);
U = zeros(size, 1);
C = zeros(size, 1);

%matrix A
for i=1:size
    A(i, i)=alpha;
end

for i=1:size-1
    A(i, i+1)=1;
end

for j=1:size-1
    A(j+1, j)=1;
end

for i=imax-2:imax-2:size-(imax-2)
    A(i, i+1)=0;
    A(i+1, i)=0;
end

for i=1:size-(imax-2)
    A(i, i+(imax-2))=beta^2;
end

for j=1:size-(imax-2)
    A(j+(imax-2), j)=beta^2;
end

%matrix C
for j=1:jmax-2
    for i=1:imax-2
        C(i+(j-1)*(imax-2))=-(u(i+2, j+1)+u(i, j+1)+beta^2*u(i+1, j+2)+beta^2*u(i+1, j));
    end
end

U = A\C; %gauss elimination

n = 1;

for j=2:jmax-1
    for i=2:imax-1
        u(i, j) = U(n);
        n=n+1;
    end
end

end