% Copyright, 2012, Michael Andersen, DIKU.
clear all;
close all;
clc;

tol_rel = 0.0;
tol_abs = 0.0;

n        = 30;  % Dimension of problem
max_iter = 100; % Maximum number of solver iterations

A      = make_contact_matrix(n);
[x0,b] = make_lcp(A,0.25);

[z e i f conv m] = ip_lcp(A, b, x0, max_iter, tol_rel, tol_abs, true);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s"\n', abs(A*z+b)'*z, e, i, m);

clf;
set(gca,'FontSize',18);
h1 = semilogy(conv, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
title('Interior point solver: contact LCP','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Complementarity measure x''*y','FontSize',18);
print('-f1','-depsc2','output/ip_solver');