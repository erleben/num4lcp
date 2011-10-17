% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;
A      = make_fluid_matrix(10);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));

tol_rel  = 0.0;
tol_abs  = 0.0;
max_iter = 30;

    
[z1 e1 i1 f1 conv1 m1] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random', true);
display( strcat( 'Parameter test: Fischer = ', m1(f1), ' in ', num2str(i1), ' iterations' ));

[z2 e2 i2 f2 conv2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
display( strcat( 'Parameter test: Fischer = ', m2(f2), ' in ', num2str(i2), ' iterations' ));

[z3 e3 i3 f3 conv3 m3] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'zero', true);
display( strcat( 'Parameter test: Fischer = ', m3(f3), ' in ', num2str(i3), ' iterations' ));

[z4 e4 i4 f4 conv4 m4] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'approximation', true);
display( strcat( 'Parameter test: Fischer = ', m4(f4), ' in ', num2str(i4), ' iterations' ));

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, ':o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-.s','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '--d','LineWidth',2,'Color',[0.7 0.1, 0.7]);
title('Comparison of Newton system solver: fluid LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2, h3, h4], 'Random','Perturbation','Zero', 'Approximation');
hold off;
print('-f1','-depsc2','output/comparison_newton_system_solver_fluid');

A      = make_contact_matrix(30);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z1 e1 i1 f1 conv1 m1] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random', true);
display( strcat( 'Parameter test: Fischer = ', m1(f1), ' in ', num2str(i1), ' iterations' ));

[z2 e2 i2 f2 conv2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
display( strcat( 'Parameter test: Fischer = ', m2(f2), ' in ', num2str(i2), ' iterations' ));

[z3 e3 i3 f3 conv3 m3] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'zero', true);
display( strcat( 'Parameter test: Fischer = ', m3(f3), ' in ', num2str(i3), ' iterations' ));

[z4 e4 i4 f4 conv4 m4] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'approximation', true);
display( strcat( 'Parameter test: Fischer = ', m4(f4), ' in ', num2str(i4), ' iterations' ));

figure(2)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, ':o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-.s','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '--d','LineWidth',2,'Color',[0.7 0.1, 0.7]);
title('Comparison of Newton system solver: contact LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2, h3, h4], 'Random','Perturbation','Zero', 'Approximation');
hold off;
print('-f2','-depsc2','output/comparison_newton_system_solver_contact');
