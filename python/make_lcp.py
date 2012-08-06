import numpy as np
# from print_vector import print_vector

def make_lcp(A,F):
    """
    Copyright 2012, Michael Andersen, DIKU, michael (at) diku (dot) dk
    """

    # [x, b] = MAKE_LCP(A): Make LCP
    #
    # INPUT:
    #
    #   A -- The coefficient matrix in the LCP
    #   F -- The fraction of zero-values in the x-solution.
    #
    # OUTPUT:
    #
    #   x -- A solution for the LCP problem.
    #   b -- The right hand side vector
    #
    # Port of Kenny Erleben, DIKU 2011 Matlab code to Python

    ##### Get of number of variable ##########
    N = np.size(A,0) # Pick a dimension, it should be NxN

    ##### Generate a random LCP solution ##########
    x = np.random.uniform(0,1,(N,1))
    x[x < F] = 0

    # x = np.real_if_close(x)

    ##### Generate a right-hand-side vector that is ##########
    ##### consistent with a random solution         ##########
    b = np.zeros((N,1))
    s = np.real_if_close(np.dot(A,x))
    b[x>0] = -s[x>0]

    return (x, b)


