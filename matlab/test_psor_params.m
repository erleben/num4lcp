% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;

lambda1 = 0.2;
lambda2 = 0.4;
lambda3 = 0.6;
lambda4 = 0.8;
lambda5 = 1.0;
lambda6 = 1.2;
lambda7 = 1.4;
lambda8 = 1.6;
lambda9 = 1.8;

max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;

A      = make_fluid_matrix(10);
[x, b] = make_lcp(A,0.25);
x0     = zeros(size(x));
    
[z1 e1 i1 f1 conv1 m1] = psor(A, b, x0, lambda1, max_iter, tol_rel, tol_abs, true);
[z2 e2 i2 f2 conv2 m2] = psor(A, b, x0, lambda2, max_iter, tol_rel, tol_abs, true);
[z3 e3 i3 f3 conv3 m3] = psor(A, b, x0, lambda3, max_iter, tol_rel, tol_abs, true);
[z4 e4 i4 f4 conv4 m4] = psor(A, b, x0, lambda4, max_iter, tol_rel, tol_abs, true);
[z5 e5 i5 f5 conv5 m5] = psor(A, b, x0, lambda5, max_iter, tol_rel, tol_abs, true);
[z6 e6 i6 f6 conv6 m6] = psor(A, b, x0, lambda6, max_iter, tol_rel, tol_abs, true);
[z7 e7 i7 f7 conv7 m7] = psor(A, b, x0, lambda7, max_iter, tol_rel, tol_abs, true);
[z8 e8 i8 f8 conv8 m8] = psor(A, b, x0, lambda8, max_iter, tol_rel, tol_abs, true);
[z9 e9 i9 f9 conv9 m9] = psor(A, b, x0, lambda9, max_iter, tol_rel, tol_abs, true);

display( strcat( 'Parameter test: final PGS state = ', m1(f1), ' in ', num2str(i1), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m2(f2), ' in ', num2str(i2), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m3(f3), ' in ', num2str(i3), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m4(f4), ' in ', num2str(i4), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m5(f5), ' in ', num2str(i5), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m6(f6), ' in ', num2str(i6), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m7(f7), ' in ', num2str(i7), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m8(f8), ' in ', num2str(i8), ' iterations' ));
display( strcat( 'Parameter test: final PGS state = ', m9(f9), ' in ', num2str(i9), ' iterations' ));

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(conv1, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(conv2, '-o','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = semilogy(conv3, '-s','LineWidth',2,'Color',[0.1 0.1, 0.7]);
h4 = semilogy(conv4, '-d','LineWidth',2,'Color',[0.7 0.1, 0.7]);
h5 = semilogy(conv5, ':','LineWidth',2,'Color',[0.1 0.7, 0.7]);
h6 = semilogy(conv6, '-.','LineWidth',2,'Color',[0.7 0.7, 0.1]);
h7 = semilogy(conv7, '--','LineWidth',2,'Color',[0.7 0.3, 0.3]);
h8 = semilogy(conv8, '-*','LineWidth',2,'Color',[0.3 0.7, 0.3]);
h9 = semilogy(conv9, '-v','LineWidth',2,'Color',[0.3 0.3, 0.7]);
title('Parameter study of PSOR','FontSize',18);
xlabel('Iterations','FontSize',18);
ylabel('Merit value','FontSize',18);
legend([h1, h2, h3, h4 h5 h6 h7 h8 h9], '\lambda = 0.2','\lambda = 0.4','\lambda = 0.6','\lambda = 0.8','\lambda = 1.0','\lambda = 1.2','\lambda = 1.4','\lambda = 1.6','\lambda = 1.8');
hold off;
print('-f1','-depsc2','output/param_psor');
