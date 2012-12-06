% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;

k1 = 1;
k2 = 2;
k3 = 3;
k4 = 4;
k5 = 5;

max_iter = 25;
tol_rel  = 0.0;
tol_abs  = 0.0;

A      = make_fluid_matrix(10);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z1 e1 i1 f1 conv1 m1] = multilevel(A, b, x0, k1, max_iter, tol_rel, tol_abs, true);
display( [ 'Parameter test: multilevel = ' m1 ' in ' num2str(i1) ' iterations' ] );

[z2 e2 i2 f2 conv2 m2] = multilevel(A, b, x0, k2, max_iter, tol_rel, tol_abs, true);
display( [ 'Parameter test: multilevel = ' m2 ' in ' num2str(i2) ' iterations' ] );

[z3 e3 i3 f3 conv3 m3] = multilevel(A, b, x0, k3, max_iter, tol_rel, tol_abs, true);
display( [ 'Parameter test: multilevel = ' m3 ' in ' num2str(i3) ' iterations' ] );

[z4 e4 i4 f4 conv4 m4] = multilevel(A, b, x0, k4, max_iter, tol_rel, tol_abs, true);
display( [ 'Parameter test: multilevel = ' m4 ' in ' num2str(i4) ' iterations'] );

[z5 e5 i5 f5 conv5 m5] = multilevel(A, b, x0, k5, max_iter, tol_rel, tol_abs, true);
display( [ 'Parameter test: multilevel = ' m5 ' in ' num2str(i5) ' iterations' ] );

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, '--o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-.d','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, ':s','LineWidth',2,'Color',[0.7 0.1, 0.7]);
h5 = semilogy(conv5, '-v','LineWidth',2,'Color',[0.1 0.7, 0.7]);
title('Parameter study of Multilevel','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2, h3, h4 h5], 'k = 1','k=2','k=3','k=4','k=5');
hold off;
print('-f1','-depsc2','output/param_multilevel');
