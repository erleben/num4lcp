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
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''random'' solver.\n', abs(A*z1+b)'*z1, e1, i1, m1);

[z2 e2 i2 f2 conv2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''perturbation'' solver.\n', abs(A*z2+b)'*z2, e2, i2, m2);

[z3 e3 i3 f3 conv3 m3] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'zero', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''zero'' solver.\n', abs(A*z3+b)'*z3, e3, i3, m3);

[z4 e4 i4 f4 conv4 m4] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'approximation', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''approximation'' solver.\n', abs(A*z4+b)'*z4, e4, i4, m4);

[z5 e5 i5 f5 conv5 m5] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'penalized', true, 0.90);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''penalized lambda = 0.90'' solver.\n', abs(A*z5+b)'*z5, e5, i5, m5);

[z6 e6 i6 f6 conv6 m6] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'penalized', true, 1.00);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''penalized lambda = 1.00'' solver.\n', abs(A*z6+b)'*z6, e6, i6, m6);

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, ':o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-.s','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '--d','LineWidth',2,'Color',[0.7 0.1, 0.7]);
h5 = semilogy(conv5, ':+','LineWidth',2,'Color',[0.1 0.7, 0.7]);
h6 = semilogy(conv6, '-p','LineWidth',2,'Color',[0.7 0.7, 0.1]);
title('Comparison of Newton system solver: fluid LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2, h3, h4 h5 h6], 'Random','Perturbation','Zero', 'Approximation', 'Penalized \lambda = 0.90', 'Penalized \lambda = 1.00');
hold off;
print('-f1','-depsc2','output/comparison_newton_system_solver_fluid');

A      = make_contact_matrix(50);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z1 e1 i1 f1 conv1 m1] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''random'' solver.\n', abs(A*z1+b)'*z1, e1, i1, m1);

[z2 e2 i2 f2 conv2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''perturbation'' solver.\n', abs(A*z2+b)'*z2, e2, i2, m2);

[z3 e3 i3 f3 conv3 m3] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'zero', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''zero'' solver.\n', abs(A*z3+b)'*z3, e3, i3, m3);

[z4 e4 i4 f4 conv4 m4] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'approximation', true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''approximation'' solver.\n', abs(A*z4+b)'*z4, e4, i4, m4);

[z5 e5 i5 f5 conv5 m5] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'penalized', true, 0.90);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''penalized lambda = 0.90'' solver.\n', abs(A*z5+b)'*z5, e5, i5, m5);

[z6 e6 i6 f6 conv6 m6] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'penalized', true, 1.00);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s using ''penalized lambda = 1.00'' solver.\n', abs(A*z6+b)'*z6, e6, i6, m6);

figure(2)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, ':o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-.s','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '--d','LineWidth',2,'Color',[0.7 0.1, 0.7]);
h5 = semilogy(conv5, ':+','LineWidth',2,'Color',[0.1 0.7, 0.7]);
h6 = semilogy(conv6, '-p','LineWidth',2,'Color',[0.7 0.7, 0.1]);

title('Comparison of Newton system solver: contact LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1 h2 h3 h4 h5 h6], 'Random','Perturbation','Zero', 'Approximation', 'Penalized \lambda = 0.90', 'Penalized \lambda = 1.00' );
hold off;
print('-f2','-depsc2','output/comparison_newton_system_solver_contact');
