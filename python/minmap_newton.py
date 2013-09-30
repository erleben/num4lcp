import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import dsolve

# from pyamg.relaxation import block_gauss_seidel

import sys
# from print_vector import print_vector

def minmap_newton(A, b, x, max_iter=0, tol_rel=0.00001, tol_abs=np.finfo(np.float64).eps*10, solver='random', profile=True):
    """
    Copyright 2012, Michael Andersen, DIKU, michael (at) diku (dot) dk
    """

    # Human readable dictionary of exit messages
    msg = {1 : 'preprocessing',  # flag = 1
           2 : 'iterating',      # flag = 2
           3 : 'relative',       # flag = 3
           4 : 'absolute',       # flag = 4
           5 : 'stagnation',     # flag = 5
           6 : 'local minima',   # flag = 6
           7 : 'nondescent',     # flag = 7
           8 : 'maxlimit',       # flag = 8
           }

    # We use N as a reference for the usage of asserts throughout the
    # program. We assume that everything is column vector of shape
    # (N,1) if not the asserts will catch them.
    N = np.size(b)
    flag = 1

    assert x.shape == (N,1), 'x0 is not a column vector, it has shape: ' + repr(x.shape)
    assert A.shape == (N,N), 'A is not a square matrix, it has shape: ' + repr(A.shape)
    assert b.shape == (N,1), 'b is not a column vector, it has shape: ' + repr(b.shape)

    if max_iter == 0:
        max_iter = np.floor(N/2.0)

    # Ensure sane values
    max_iter = max(max_iter,1)
    # Rest of the value should be sane

    ##### Magic constants #####
    h      = 1e-7
    alpha  = 0.5
    beta   = 0.001
    gamma  = 1e-28

    eps    = np.finfo(np.float64).eps
    rho    = np.finfo(np.float64).eps
    gmres_tol = 10*eps

    # ##### Warm start of FN using Blocked Gauss-Seidel from PyAMG #####
    # warm_start = False
    # max_warm_iterate = 5
    # max_bgs_iterate = 5
    # if warm_start:
    #     x = fischer_warm_start(A, x, b, max_warm_iterate, max_bgs_iterate)

    ##### Values needed file iterating #####
    convergence = np.zeros(max_iter+1)

    # Should use np.infty.
    err     = 1e20
    iterate = 1
    flag    = 2

    while iterate <= max_iter:
        y = np.dot(A,x) + b
        assert y.shape == (N,1), 'y is not a column vector, it has shape: ' + repr(y.shape)
        assert np.all(np.isreal(y)), 'y is not real'
        # Calculate the minimum map column vector.
        H = minmap(x,y)
        assert H.shape == (N,1), 'H is not a column vector, it has shape: ' + repr(H.shape)
        assert np.all(np.isreal(H)), 'H is not real'
        old_err = err
        # Calculate merit value, error
        err = 0.5*np.dot(H.T,H)
        assert err.shape == (1,1), 'err is not a scalar, it has shape: ' + repr(err.shape)
        assert np.isreal(err), 'err is not real'

        if profile:
            convergence[iterate-1] = err

        ##### Test the stopping criterias used #####
        rel_err = np.abs(err-old_err)/np.abs(old_err)

        if rel_err < tol_rel:
            flag = 3
            break

        if err < tol_abs:
            flag = 4
            break

        ##### Solving the Newton system
        restart = min(N, 20) # Number of iterates done before Restart
                             # for GMRES should restart
        S = np.where(y < x)
        J = np.identity(N)
        J[S,:] = A[S,:]
        dx = np.zeros((N,1))
        dx = gmres(J, (-H), tol=gmres_tol, restart=restart)[0].reshape(N,1)

        assert dx.shape == (N,1), 'dx is not a column vector, it has shape: ' + repr(dx.shape)
        assert np.all(np.isreal(dx)), 'dx is not real'

        nabla_H = np.dot(H.T, J)
        # Ensure nabla_H is a column vector
        nabla_H = nabla_H.reshape(N,1)
        assert nabla_H.shape == (N,1), 'nabla_H is not a column vector, it has shape: ' + repr(nabla_H.shape)
        assert np.all(np.isreal(nabla_H)), 'nabla_H is not real'

        # Tests whether the search direction is below machine
        # precision.
        if np.max(np.abs(dx)) < eps:
            flag = 5
            print "*** Search direction below machine precision, choosing gradient as search direction."
            dx = -nabla_H

        # Test whether we are stuck in a local minima
        if np.linalg.norm(nabla_H) < tol_abs:
            flag = 6
            break

        # Test whether our direction is a sufficient descent direction
        if np.dot(nabla_H.T,dx) > -rho*(np.dot(dx.T, dx)):
            # Otherwise we should try gradient direction instead.
            print "*** Non descend direction at iterate " + repr(iterate) + ", choosing gradient as search direction."
            dx = -nabla_H

        ##### Armijo backtracking combined with a projected line-search #####
        tau = 1.0
        f_0 = err
        grad_f = beta*np.dot(nabla_H.T,dx)

        x_k = x[:]
        assert x_k.shape == (N,1), 'x_k is not a column vector, it has shape: ' + repr(x_k)
        assert np.all(np.isreal(x_k)), 'x_k is not real'
    
        # Perform backtracking line search
        while True:
            x_k = np.maximum(0, x + dx*tau)
            assert x_k.shape == (N,1), 'x_k is not a column vector, it has shape: ' + repr(x_k.shape)
            assert np.all(np.isreal(x_k)), 'x_k is not real'
            y_k = np.dot(A,x_k)+b
            assert y_k.shape == (N,1), 'y_k is not a column vector, it has shape: ' + repr(y_k.shape)
            assert np.all(np.isreal(y_k)), 'y_k is not real'
            H_k = minmap(y_k,x_k)
            assert H_k.shape == (N,1), 'H_k is not a column vector, it has shape: ' + repr(H_k.shape)
            assert np.all(np.isreal(H_k)), 'H_k is not real'
            f_k = 0.5*(np.dot(H_k.T,H_k))
            # Test Armijo condition for sufficient decrease
            if f_k <= f_0 + tau*grad_f:
                break
            # Test whether the stepsize has become too small
            if tau*tau < gamma:
                break
            tau *= alpha

        # Update iterate with result from line search.
        x = x_k
        assert x.shape == (N,1), 'x is not a column vector, it has shape: ' + repr(x.shape)
        assert np.all(np.isreal(x)), 'x is not real.'

        # Increment iterate
        iterate += 1

    if iterate >= max_iter:
        iterate -= 1
        flag = 8

    return (x, err, iterate, flag, convergence[:iterate], msg[flag])

def minmap(x,y):
    return np.minimum(x,y)
