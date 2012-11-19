%%% Script for testing penalized Fischer-Burmeister %%%

%%% Static test
M = [0  1 0;
     0  0 1;
     0 -1 1];
 
q = [0; 0; 1];
 
x0 = ones(size(M,1), 1);

max_iter = 100;
tol_rel = 0.0001;
tol_abs = eps*10;
 
[z1 e1 i1 f1 c1 m1] = penalized_fischer_burmeister(M, q, x0, max_iter);
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s\n', abs(M*z1+q)'*z1, e1, i1, m1);
    
[z2 e2 i2 f2 c2 m2] = fischer_newton(M, q, x0, max_iter, tol_rel, tol_abs,'perturbation', true );
fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s"\n', abs(M*z2+q)'*z2, e2, i2, m2);

%%% Test with synthetic lcp problems %%%

num_test = 100;

% for i = 1:num_test
    A      = make_contact_matrix(5);
    [x, b] = make_lcp(A,0.25);
    x0     = zeros(size(x));

    [z1 e1 i1 f1 c1 m1] = penalized_fischer_burmeister(A, b, x0, max_iter, tol_rel, tol_abs, true);
    fprintf('Penalized Fischer-Burmeister:\n')
    fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s\n', abs(A*z1+b)'*z1, e1, i1, m1);
    
    [z2 e2 i2 f2 c2 m2] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', true);
    fprintf('Fischer Newton:\n')
    fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s"\n', abs(A*z2+b)'*z2, e2, i2, m2);        
% end
