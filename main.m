clear
close all

%input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=21;
jmax=41;
L=1;
H=2;

method = 1;
maxiteration = 5000;
tolerance = 0.0001;
exactsolution = 1; %exact solution을 같이 표시하려면 1, 아니면 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method
%1.FTCS  2.Richardson  3.DuFort-Frankel 4.Laasonen  5.Crank-Nicolson
%6.Beta formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%initialization%%%%%
u = initialization(imax, jmax);

dx = L/(imax-1);
dy = H/(jmax-1);

%%%%%method%%%%%
if method == 1
    u = Five_point(u, dx, dy, imax, jmax);
elseif method == 2
    [u, residual] = Jacobi(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 3
    u = D_F(u, alpha, dx, dt, imax, nmax);
elseif method == 4
    u = Laasonen(u, alpha, dx, dt, imax, nmax);
elseif method == 5
    u = C_N(u, alpha, dx, dt, imax, nmax);
elseif method == 6
    u = Beta(beta, u, alpha, dx, dt, imax, nmax);
end

contourf((0:imax-1)*dx, (0:jmax-1)*dy, u.')