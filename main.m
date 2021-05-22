clear
close all

%input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=21;
jmax=41;
L=1;
H=2;

method = 2;
maxiteration = 5000;
tolerance = 0.0001;
exactsolution = 1; %exact solution을 같이 표시하려면 1, 아니면 0
level = 10:10:90;
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

%%%%%exact solution%%%%%
if exactsolution
    uex = exact(H, L, dx, dy, imax, jmax);
end

%%%%%method%%%%%
if method == 1
    u = Five_point(u, dx, dy, imax, jmax);
elseif method == 2
    [u, residual] = Jacobi(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 3
    u = exact(H, L, dx, dy, imax, jmax);
elseif method == 4
    u = Laasonen(u, alpha, dx, dt, imax, nmax);
elseif method == 5
    u = C_N(u, alpha, dx, dt, imax, nmax);
elseif method == 6
    u = Beta(beta, u, alpha, dx, dt, imax, nmax);
end

hold on

if exactsolution
contour((0:imax-1)*dx, (0:jmax-1)*dy, uex.', level, ':', 'color', [0.5 0.5 0.5], 'linewidth', 5)
end

[C, h]=contour((0:imax-1)*dx, (0:jmax-1)*dy, u.', level, 'k', 'linewidth', 1.5);
clabel(C, h, 'backgroundcolor', 'w')

set(gcf, 'position', [100 100 300 600])