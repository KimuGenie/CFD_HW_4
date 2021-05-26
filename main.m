%main program

clear
close all

%input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=21;
jmax=41;
L=1;
H=2;

method = 5;
maxiteration = 1000;
tolerance = 0.0001;
omega = 1.7; %relaxation parameter

exactsolution = 1; %exact solution을 같이 표시하려면 1, 아니면 0

level = 10:10:90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method
%1.Five Point  2.Jacobi  3.Point Gauss-Seidel 4.Line Gauss-Seidel 
%5.Point Successive Over-Relaxation
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
    [u, residual] = PGS(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 4
    [u, residual] = LGS(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 5
    [u, residual] = PSOR(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
end

hold on

figure(1)
if exactsolution
contour((0:imax-1)*dx, (0:jmax-1)*dy, uex.', level, ':', 'color', [0.5 0.5 0.5], 'linewidth', 5)
end

[C, h]=contour((0:imax-1)*dx, (0:jmax-1)*dy, u.', level, 'k', 'linewidth', 1.5);
clabel(C, h, 'backgroundcolor', 'w')
set(gcf, 'position', [100 100 300 600])

figure(2)
plot(residual, 'k', 'linewidth', 1.5)