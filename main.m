%main program

clear
close all

%input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax=21;
jmax=41;
L=1;
H=2;

maxiteration = 5000;
tolerance = 1e-4;
omega = 1.3; %relaxation parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode = 1;
%1.contour그리기  2.method 수렴성 평가  3.method 실행시간 평가
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 4; %contour를 그릴 method
%method
%1.Five Point  2.Jacobi  3.Point Gauss-Seidel  4.Line Gauss-Seidel 
%5.Point Successive Over-Relaxation  6.Line Successive Over-Relaxation
%7.Alternating Direction Implicit method
%8.Alternating Direction Implicit Over-Relaxation method
%9.Line Jacobi

exactsolution = 1; %exact solution을 같이 표시하려면 1, 아니면 0

level = 10:10:90; %contour의 level

test = 10; %실행시간 측정 테스트 횟수
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%initialization%%%%%
u = initialization(imax, jmax);

dx = L/(imax-1);
dy = H/(jmax-1);

%%%%%exact solution%%%%%
uex = exact(H, L, dx, dy, imax, jmax);


%%%%%%%%%%%%%%%%%%%%contour 그리기%%%%%%%%%%%%%%%%%%%%
if mode == 1

residual=[];

%%%%%method%%%%%
if method == 1
    u = Five_point(u, dx, dy, imax, jmax);
elseif method == 2
    [u, residual] = Jacobi(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 3
    [u, residual] = PGS(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 4
    [u, residual] = LGSx(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 5
    [u, residual] = PSOR(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
elseif method == 6
    [u, residual] = LSORy(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
elseif method == 7
    [u, residual] = ADIxy(u, dx, dy, imax, jmax, maxiteration, tolerance);
elseif method == 8
    [u, residual] = ADIORxy(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
elseif method == 9
    [u, residual] = LJacobiy(u, dx, dy, imax, jmax, maxiteration, tolerance);
end

hold on

figure(1)
if exactsolution
contour((0:imax-1)*dx, (0:jmax-1)*dy, uex.', level, ':', 'color', [0.5 0.5 0.5], 'linewidth', 5)
end

[C, h]=contour((0:imax-1)*dx, (0:jmax-1)*dy, u.', level, 'k', 'linewidth', 1.5);
clabel(C, h, 'backgroundcolor', 'w')
set(gcf, 'position', [100 100 300 600])

disp(['iteration=', num2str(length(residual))])


%%%%%%%%%%%%%%%%%%%%수렴성 평가%%%%%%%%%%%%%%%%%%%%
elseif mode == 2
    
%%%%%method%%%%%
u0 = Five_point(u, dx, dy, imax, jmax);
[u1, residual1] = Jacobi(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u2, residual2] = PGS(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u3, residual3] = LGSx(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u4, residual4] = LGSy(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u5, residual5] = PSOR(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
[u6, residual6] = LSORx(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
[u7, residual7] = LSORy(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
[u8, residual8] = ADIxy(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u9, residual9] = ADIyx(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u10, residual10] = ADIORxy(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
[u11, residual11] = ADIORyx(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
[u12, residual12] = LJacobix(u, dx, dy, imax, jmax, maxiteration, tolerance);
[u13, residual13] = LJacobiy(u, dx, dy, imax, jmax, maxiteration, tolerance);

lin = ["k-o", "k-s", "k-^", "k-x", "k-p", "k-v", "k-*", "k-+", "k-<", "k-d"];
n=1;
ms=30;


figure(1)
semilogy(residual1, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual1), 'DisplayName', 'Jacobi'); n=n+1;
hold on
semilogy(residual2, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual2), 'DisplayName', 'PGauss'); n=n+1;
semilogy(residual3, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual3), 'DisplayName', 'LGaussx'); n=n+1;
% hold on
% semilogy(residual4, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual4), 'DisplayName', 'LGaussy'); n=n+1;
semilogy(residual5, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual5), 'DisplayName', 'PSOR'); n=n+1;
semilogy(residual6, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual6), 'DisplayName', 'LSORx'); n=n+1;
% hold on
% semilogy(residual7, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual7), 'DisplayName', 'LSORy'); n=n+1;
semilogy(residual8, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual8), 'DisplayName', 'ADIxy'); n=n+1;
% hold on
% semilogy(residual9, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual9), 'DisplayName', 'ADIyx'); n=n+1;
semilogy(residual10, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual10), 'DisplayName', 'ADIORxy'); n=n+1;
hold on
% semilogy(residual11, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual11), 'DisplayName', 'ADIORyx'); n=n+1;
% semilogy(residual12, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual12), 'DisplayName', 'LJacobix'); n=n+1;
% hold on
semilogy(residual13, lin(n), 'linewidth', 1.5, 'MarkerIndices', 1:ms:length(residual13), 'DisplayName', 'LJacobiy'); n=n+1;


legend

title("ω="+num2str(omega)+" tolerance="+num2str(tolerance), 'fontsize', 13)
xlabel 'Iteration'
ylabel 'Total change in temperature(°R)'
set(gcf, 'position', [100 300 400 300])
ylim([tolerance 10])

figure(2)
err(1)=mean(abs(u0-uex), 'all');
err(2)=mean(abs(u1-uex), 'all');
err(3)=mean(abs(u13-uex), 'all');
err(4)=mean(abs(u2-uex), 'all');
err(5)=mean(abs(u3-uex), 'all');
err(6)=mean(abs(u8-uex), 'all');
err(7)=mean(abs(u5-uex), 'all');
err(8)=mean(abs(u6-uex), 'all');
err(9)=mean(abs(u10-uex), 'all');

methodname = categorical(["Five-point", "Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);
methodname = categorical(methodname, ["Five-point", "Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);

bar(methodname, err)

title("ω="+num2str(omega)+" tolerance="+num2str(tolerance), 'fontsize', 13)
ylabel 'Mean of error(°R)'
set(gcf, 'position', [500 300 400 300])

figure(3)
iteration(1)=length(residual1);
iteration(2)=length(residual13);
iteration(3)=length(residual2);
iteration(4)=length(residual3);
iteration(5)=length(residual8);
iteration(6)=length(residual5);
iteration(7)=length(residual6);
iteration(8)=length(residual10);

methodname = categorical(["Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);
methodname = categorical(methodname, ["Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);

bar(methodname, iteration)

title("ω="+num2str(omega)+" tolerance="+num2str(tolerance), 'fontsize', 13)
ylabel 'Number of iterations'
set(gcf, 'position', [900 300 400 300])

%%%%%%%%%%%%%%%%%%%%실행시간 평가%%%%%%%%%%%%%%%%%%%%
elseif mode == 3
    
%%%%%methods%%%%%
td = zeros(test, 1);
eu = zeros(8, 1);
el = zeros(8, 1);
tm = zeros(8, 1);

methodname = categorical(["Five-point", "Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);
methodname = categorical(methodname, ["Five-point", "Jacobi", "LJacobiy", "PGauss", "LGaussx", "ADIxy", "PSOR", "LSORx", "ADIORxy"]);

for n=1:test
    tic
    ud = Five_point(u, dx, dy, imax, jmax);
    td(n)=toc;
end
tm(1) = mean(td);
eu(1) = max(td)-tm(1);
el(1) = tm(1)-min(td);

itermethods=["Jacobi", "LJacobiy", "PGS", "LGSx", "ADIxy"];
for m=1:5
    for n=1:test
    f = str2func(itermethods(m));
    tic
    [ud, rd] = f(u, dx, dy, imax, jmax, maxiteration, tolerance);
    td(n)=toc;
    end
    tm(m+1) = mean(td);
    eu(m+1) = max(td)-tm(m+1);
    el(m+1) = tm(m+1)-min(td);
end

iterormethods=["PSOR", "LSORx", "ADIORxy"];
for m=1:3
    f = str2func(iterormethods(m));
    for n=1:test
    tic
    [ud, rd] = f(u, dx, dy, imax, jmax, maxiteration, tolerance, omega);
    td(n)=toc;
    end
    tm(m+6) = mean(td);
    eu(m+6) = max(td)-tm(m+6);
    el(m+6) = tm(m+6)-min(td);
end

tm = tm.*1000; %ms로 환산
el = el.*1000;
eu = eu.*1000;

bar(methodname, tm)

hold on

er = errorbar(methodname, tm, el, eu);    
er.Color = [0 0 0];
er.LineStyle = 'none';

ylabel 'Execute time(ms)'
set(gcf, 'position', [500 300 400 300])

end