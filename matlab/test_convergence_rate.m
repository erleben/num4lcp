% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;

max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;

A      = make_fluid_matrix(10);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z1 e1 i1 f1 conv1 m1] = multilevel(A, b, x0, 2, max_iter, tol_rel, tol_abs, true);
display( strcat( 'Convergence test: Multilevel = ',     m1(f1), ' in ', num2str(i1), ' iterations'));

[z2 e2 i2 f2 conv2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
display( strcat( 'Convergence test: Fischer = ', m2(f2), ' in ', num2str(i2), ' iterations'));

[z3 e3 i3 f3 conv3 m3] = minmap_newton(A, b, x0, max_iter, tol_rel, tol_abs, true);
display( strcat( 'Convergence test: Min map = ', m3(f3), ' in ', num2str(i3), ' iterations'));

[z4 e4 i4 f4 conv4 m4] = pgs(A, b, x0, max_iter, tol_rel, tol_abs, true);
display( strcat( 'Convergence test: PGS = ',            m4(f4), ' in ', num2str(i4), ' iterations'));

[z5 e5 i5 f5 conv5 m5] = psor(A, b, x0, 1.4, max_iter, tol_rel, tol_abs, true);
display( strcat( 'Convergence test: PSOR = ',            m5(f5), ' in ', num2str(i5), ' iterations'));

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
title('Multilevel Convergence Rate: fluid LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
hold off;
print('-f1','-depsc2','output/convergence_multilevel_fluid');

figure(2)
clf;
set(gca,'FontSize',18);
h2 = semilogy(conv2, ':','LineWidth',2,'Color',[0.1 0.7, 0.1]);
grid on;
hold on;
h3 = semilogy(conv3, '-.','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '--','LineWidth',2,'Color',[0.7 0.1, 0.7]);
h5 = semilogy(conv5, ':','LineWidth',2,'Color',[0.1 0.7, 0.7]);
title('Convergence Rate: fluid LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h2, h3, h4 h5], 'Fischer','Min Map', 'PGS', 'PSOR');
hold off;
print('-f2','-depsc2','output/convergence_fluid');

A      = make_contact_matrix(50);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z6 e6 i6 f6 conv6 m6] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
display( strcat( 'Convergence test: Fischer = ', m6(f6), ' in ', num2str(i1), ' iterations') );

[z7 e7 i7 f7 conv7 m7] = minmap_newton(A, b, x0, max_iter, tol_rel, tol_abs, true);
display( strcat( 'Convergence test: Min map = ', m7(f7), ' in ', num2str(i2), ' iterations') );

figure(3)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv6, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv7, '-','LineWidth',2,'Color',[0.1 0.7, 0.1]);
title('Convergence Rate: contact LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2],'Fischer','Min Map');
hold off;
print('-f3','-depsc2','output/convergence_contact');
