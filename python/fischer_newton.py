import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator
import sys
from print_vector import print_vector

def fischer_newton(A, b, x, max_iter=0, tol_rel=0.00001, tol_abs=np.finfo(np.float64).eps*10, solver='random', profile=True):
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

    ##### Values needed file iterating #####
    convergence = np.zeros(max_iter+1)
    # Should use np.infty.
    err     = 1e20
    iterate = 1
    flag    = 2

    grad_steps = 0

    while iterate <= max_iter:

        y = np.dot(A,x)+b
        assert y.shape == (N,1), 'y is not a column vector, it has shape: ' + repr(y.shape)
        assert np.all(np.isreal(y)), 'y is not real'

        # Calculate the fischer function value
        phi = fischer(y,x)
        assert phi.shape == (N,1), 'phi is not a column vector, it has shape: ' + repr(phi.shape)
        assert np.all(np.isreal(phi)), 'phi is not real'
        old_err = err
        # Calculate merit value, error
        err = 0.5*np.dot(phi.T,phi)
        assert err.shape == (1,1), 'err is not a scalar, it has shape: ' + repr(err.shape)
        assert np.isreal(err), 'err is not real'
        
        if profile:
            convergence[iterate-1] = err

        ##### Test the stopping criterias used #####
        if np.abs(err-old_err)/np.abs(old_err) < tol_rel:
            flag = 3
            break

        if err < tol_abs:
            flag = 4
            break

        ##### Solving the Newton system
        dx = np.zeros((N,1))

        # Bitmask of True/False values where the conditions are
        # satisfied, that is where phi and x are near
        # singular/singular.
        S = np.logical_and(np.abs(phi) < gamma, np.abs(x) < gamma)
        # Bitmask of nonsingular values
        I = np.logical_not(S)
        idx = np.where(I)

        restart = min(np.size(A,0), 10) # Number of iterates done
                                        # before Restart for GMRES
                                        # should restart
        
        # if 'random' == str.lower(solver): # Works on full system
        #     q = np.random.rand(N,1) / np.sqrt(2) - 1.0
        #     p = np.random.rand(N,1) / np.sqrt(2) - 1.0
        #     # J = np.zeros((N,N))
        #     J = sps.lil_matrix((N,N))
        #     # Reinitialize dx, to sparse matrix, otherwise we get a
        #     # broadcast error when assigning to dx[idx]
        #     dx = sps.lil_matrix((N,1))
            
        #     q[idx] = (y[idx]/((y[idx]**2+x[idx]**2)**0.5))-1.0
        #     p[idx] = (x[idx]/((y[idx]**2+x[idx]**2)**0.5))-1.0
            
        #     # J needs to be sparse here, otherwise the following
        #     # assignment does not change J for some weird reason.
        #     J[idx[0],idx[0]] = np.diag(p[idx])*np.identity(np.shape(A[idx])[0]) + np.dot(np.diag(q[idx]),A[idx[0],:][:,idx[0]])
        #     # Convert it into 'dense' array form, to enable assert on shape
        #     J = J.toarray()
        #     assert J.shape == (N,N), 'J is not a square matrix, it has shape: ' + repr(J.shape)
        #     assert np.all(np.isreal(J)), 'J is not real'

        #     dx[idx] = gmres(J[idx[0],:][:,idx[0]], (-phi[idx]), tol=1e-20, restart=restart)[0].reshape(idx[0].size,1)
        #     dx = dx.toarray()
        #     dx = dx.reshape(-1,1)
        #     assert dx.shape == (N,1), 'dx is not a column vector, it has shape: ' + repr(dx.shape)
        #     assert np.all(np.isreal(dx)), 'dx is not real'

        if 'perturbation' == str.lower(solver): # Works on full system

            px = np.zeros(x.shape)
            px[:] = x # Copy x
            assert px.shape == (N,1), 'px is not a column vector, it has shape: ' + repr(px.shape)
            assert np.all( np.isreal(px) ), 'px is not real'
            direction = np.zeros(x.shape)
            direction[:] = np.sign(x)
            direction[ direction == 0 ] = 1
            px[S] = gamma * direction[S]
            p = np.zeros(px.shape)
            p = (px / (np.sqrt(y**2+px**2)))-1
            assert p.shape == (N,1), 'p is not a column vector, it has shape: ' + repr(p.shape)
            assert np.all(np.isreal(p)), 'p is not real'

            # Convert np.nan in p to 0.0
            # p = np.nan_to_num(p)

            q = np.zeros(y.shape)
            q = (y/(np.sqrt(y**2+px**2)))-1
            assert q.shape == (N,1), 'q is not a column vector, it has shape: ' + repr(q.shape)
            assert np.all(np.isreal(q)), 'q is not reeal'
            # Convert np.nan in q to 0.0
            # q = np.nan_to_num(q)

            J = np.dot(np.diag(p[:,0]),np.identity(N)) + np.dot(np.diag(q[:,0]),A)
            assert J.shape == (N,N), 'J is not a square matrix, it has shape: ' + repr(J.shape)
            assert np.all(np.isreal(J)), 'J is not real'

            # Reciprocal conditioning number of J
            rcond = 1/np.linalg.cond(J)
            if np.isnan(rcond):
                rcond = np.nan_to_num(rcond)

            # Check conditioning of J, if bad use pinv to solve the system.
            if rcond < 2*eps:
                dx = np.dot(np.linalg.pinv(J), (-phi))
            else:
                dx = np.linalg.solve(J,(-phi))

            assert dx.shape == (N,1), 'dx is not a column vector, it has shape: ' + repr(dx.shape)
            assert np.all(np.isreal(dx)), 'dx is not real'

        elif 'zero' == str.lower(solver):

            J = sps.lil_matrix((N,N))
            dx = sps.lil_matrix((N,1))

            q = (y[I]/((y[I]**2+x[I]**2)**(0.5)))-1
            p = (x[I]/((y[I]**2+x[I]**2)**(0.5)))-1

            J[idx[0],idx[0]] = np.diag(p)*np.identity(A[idx].shape[0]) + np.dot(np.diag(q),A[idx[0],:][:,idx[0]])
            J = J.toarray()

            # Call GMRES (This is migthy slow, might consider another
            # approach)
            dx[idx] = gmres(J[idx[0],:][:,idx[0]], (-phi[idx]), tol=1e-20, restart=restart)[0].reshape(idx[0].size,1)

            dx = dx.toarray()
            dx = dx.reshape(-1,1)

            assert dx.shape == (N,1), 'dx is not a column vector, it has shape: ' + repr(dx.shape)

        # elif 'approximation' == str.lower(solver):
            
        #     # J = np.zeros((N,N))
        #     J = sps.lil_matrix((N,N))
        #     q = y[idx]/((y[idx]**2+x[idx]**2)**(0.5))-1
        #     p = x[idx]/((y[idx]**2+x[idx]**2)**(0.5))-1

        #     J[idx[0],idx[0]] = np.diag(p)*np.identity(A[idx].shape[0]) + np.dot(np.diag(q[:,0]),A[idx[0],:][:,idx[0]])

        #     # Call GMRES, how do we handle a function handle for gmres?
            
        #     fun = lambda dx: (fischer(np.dot(A[idx[0],:][:,idx[0]],(x[idx]+dx[idx]*h)) + b[idx],x[idx] + dx*h) - phi[idx]) / h

        #     LO = LinearOperator((N,1), matvec = fun)

        #     dx[idx] = gmres(LO, )

        else:
            sys.stderr.write('Unknown solver method for Newton subsystem')

        nabla_phi = np.dot(phi.T, J)
        # We could use nabla_phi = nabla_phi.reshape(nabla_phi.size,
        # 1) instead this creates a column vector no matter what.
        nabla_phi = nabla_phi.T
        assert nabla_phi.shape == (N,1), 'nabla_phi is not a column vector, it has shape: ' + repr(nabla_phi.shape)
        assert np.all(np.isreal(nabla_phi)), 'nabla_phi is not real'
        # Tests whether the search direction is below machine
        # precision.
        if np.max(np.abs(dx)) < eps:
            flag = 5
            # We could try gradient direction instead.
            # dx = nabla_phi
            break

        # Test whether we are stuck in a local minima
        if np.linalg.norm(nabla_phi) < tol_abs:
            flag = 6
            break

        # Test whether our direction is a sufficient descent direction
        if np.dot(nabla_phi.T,dx) > -rho*(np.dot(dx.T, dx)):
            flag = 7
            # Otherwise we should try gradient direction instead.
            # dx = nabla_phi
            break

        ##### Armijo backtracking combined with a projected line-search #####
        tau = 1.0
        f_0 = err
        grad_f = beta*np.dot(nabla_phi.T,dx)

        x_k = x[:]
        assert x_k.shape == (N,1), 'x_k is not a column vector, it has shape: ' + repr(x_k)
        assert np.all(np.isreal(x_k)), 'x_k is not real'
    
        # Perform line search
        while True:

            x_k = np.maximum(0, x + dx*tau)
            assert x_k.shape == (N,1), 'x_k is not a column vector, it has shape: ' + repr(x_k.shape)
            assert np.all(np.isreal(x_k)), 'x_k is not real'

            y_k = np.dot(A,x_k)+b
            assert y_k.shape == (N,1), 'y_k is not a column vector, it has shape: ' + repr(y_k.shape)
            assert np.all(np.isreal(y_k)), 'y_k is not real'

            phi_k = fischer(y_k,x_k)
            assert phi_k.shape == (N,1), 'phi_k is not a column vector, it has shape: ' + repr(phi_k.shape)
            assert np.all(np.isreal(phi_k)), 'phi_k is not real'

            f_k = 0.5*(np.dot(phi_k.T,phi_k))

            # Test Armijo condition for sufficient decrease
            if np.all(f_k <= f_0 + tau*grad_f):
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
        iterate +=1

    if iterate >= max_iter:
        iterate -= 1
        flag = 8

    return (x, err, iterate, flag, convergence[:iterate], msg[flag])

def fischer(y,x):
    return np.sqrt(y**2+x**2) - y - x

# def fischerLinearOperator(A,x,dx,phi,h,idx):
#     fischer(np.dot(A[idx[0],:][:,idx[0]],(x[idx]+dx[idx]*h)),x[idx]+dx*h)-phi[idx]
