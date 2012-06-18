import numpy as np
from make_blocked_matrix import make_blocked_matrix

# Copyright 2012, Michael Andersen, DIKU

def make_contact_matrix(K):
    # [A] = MAKE_CONTACT_MATRIX(K): Generates a contact matrix that has
    # the same fill-in properties as that of the Stewart-Trinkle contact
    # formulation using a 2D friction cone.
    #
    # INPUT:
    #
    #   K -- Number of contacts.
    #
    # OUTPUT:
    #
    #   A -- The contact coefficient matrix.
    #
    # TODO LIST:
    #   
    #
    # Ported Matlab code, originally made by Kenny Erleben, DIKU 2011

    ##### Generate a random block matrix ##########
    B = 6
    F = 0.5
    A = make_blocked_matrix(B,K,F)

    for bi in np.arange(0,K):
        for bj in np.arange(0,K):
            i = bi*B
            j = bj*B
            if not A[i,j] == 0:
                A[i+5,j] = 0.5
                A[i+5,j+5] = 0
                A[i, j+4] = 0
                A[i+5,j+1:j+4] = -1
                A[i+1:i+4,j+5] = 1

    return A
