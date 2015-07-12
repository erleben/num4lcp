


# tol_rel = 0.0;
# tol_abs = 0.0;

# n        = 30;  % Dimension of problem
# max_iter = 100; % Maximum number of solver iterations

# A      = make_contact_matrix(n);
# [x0,b] = make_lcp(A,0.25);

# [z e i f conv m] = ip_lcp(A, b, x0, max_iter, tol_rel, tol_abs, true);
# fprintf('Accuracy = %2.5e\tErr = %2.5e iter: %d msg: %s"\n', abs(A*z+b)'*z, e, i, m);

# clf;
# set(gca,'FontSize',18);
# h1 = semilogy(conv, '-x','LineWidth',2,'Color',[0.7 0.1, 0.1]);
# grid on;
# title('Interior point solver: contact LCP','FontSize',18);
# xlabel('Iterations','FontSize',18);
# ylabel('Complementarity measure x''*y','FontSize',18);
# print('-f1','-depsc2','output/ip_solver');

import numpy as np
from scipy.io import savemat, loadmat
from ip_lcp import ip_lcp
from make_lcp import make_lcp
from make_contact_matrix import make_contact_matrix

import time

tol_rel = 0.0
tol_abs = 0.0
n = 30        # Problem size
max_iter = 30 # Maximum number of iterations

A = make_contact_matrix(n)
x, b = make_lcp(A,0.25)
x0 = np.zeros(np.shape(x))

x1, err1, iter1, flag1, conv1, msg1 = ip_lcp(A, b, x0, max_iter, tol_rel, tol_abs, 'central_trajectory', True)

print "> central trajectory: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x1)+b)),x1)), err1, iter1, msg1)

# x2, err2, iter2, flag2, conv2, msg2 = ip_lcp(A, b, x0, max_iter, tol_rel, tol_abs, 'potential_reduction', True)

# print "> potential reduction: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x2)+b)),x2)), err2, iter2, msg2)
