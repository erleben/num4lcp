import numpy as np
from scipy.io import savemat, loadmat
from fischer_newton import fischer_newton
from make_lcp import make_lcp
from make_contact_matrix import make_contact_matrix

import time

tol_rel = 0.0
tol_abs = 0.0
problem_size = 30
max_iter = 30

A = make_contact_matrix(problem_size)
x, b = make_lcp(A,0.25)
x0 = np.zeros(np.shape(x))

x1, err1, iter1, flag1, conv1, msg1 = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random', True)
print "> Random: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x1)+b)),x1)), err1, iter1, msg1)

savemat('conv_random_%d' % problem_size, mdict={'convergence_random' : conv1})

x2, err2, iter2, flag2, conv2, msg2 = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', True)
print "> Perturbation: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x2)+b)),x2)), err2, iter2, msg2)

savemat('conv_perturbation_%d' % problem_size, mdict={'convergence_random' : conv2})

x3, err3, iter3, flag3, conv3, msg3 = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'zero', True)
print "> Zero: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x3)+b)),x3)), err3, iter3, msg3)

savemat('conv_zero_%d' % problem_size, mdict={'convergence_random' : conv3})

# x4, err4, iter4, flag4, conv4, msg4 = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'approximation', True)
# print "> Zero: Accuracy = %2.5e Err = %2.5e iter: %d msg: %s" % (np.abs(np.dot(np.transpose((np.dot(A,x4)+b)),x4)), err4, iter4, msg4)
