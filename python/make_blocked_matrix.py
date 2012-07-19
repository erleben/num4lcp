import numpy as np

def make_blocked_matrix(B,K,F):
    """
    Copyright 2012, Michael Andersen, DIKU, michael (at) diku (dot) dk
    """

    # A = MAKE_BLOCKED_MATRIX(B,K,F): Make a random blocked matrix.
    #
    # INPUT:
    #
    #   B -- Block size
    #   K -- Number of blocks
    #   F -- Fillin protentage (number between 0 and 1)
    #
    # OUTPUT:
    #
    #   A -- A blocked symmetric positive definite matrix
    #
    # Port of Kenny Erleben, DIKU 2011 Matlab code into Python code

    N = B*K
    # Better safe than sorry.
    R = np.random.uniform(0,1,(N,N))
    A = 0.5 * (R + R.T)

    for i in np.arange(0,K):
        i_from   = i*B
        i_center = i_from + np.floor(B/2.0)
        i_to     = i_from + B

        for j in np.arange(0,K):
            j_from   = j*B
            j_center = j_from + np.floor(B/2.0)
            j_to     = j_from + B

            if A[i_center, j_center] > F:
                A[i_from:i_to,j_from:j_to] = np.zeros((B,B))
        
    w, _ = np.linalg.eig(A)
    offset = np.min(w)

    if offset <= 0:
        A = A + np.identity(N)*(0.5-offset)

    return A

