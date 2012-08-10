import numpy as np
import scipy.sparse as sps
import scipy.linalg as spl
import scipy.sparse.linalg as spsl

from ArmijoLineSearch import ArmijoLineSearch

"""
Copyright 2012, Michael Andersen, DIKU. michael (at) diku (dot) dk
"""

# We might want to delay the imports untill we actually need them.
from pyamg.relaxation import block_gauss_seidel

class NonSquareException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)
        
class InvalidShapeException(NonSquareException):
    pass

class FischerNewtonSolver:
    """
    Fischer Newton's method for solving Linear Complementarity problems

    Usage:
      fns = FischerNewtonSolver()

      (x, err, iterate, flag, convergence, msg) = fns.solve(A, b, x0, max_iter=None, tol_rel=0.0001, tol_abs=10*np.finfo(np.float64).eps, subsolver=None, profile=True, warmstart=None, gradient=None, searchmethod=None)

      @params
        A             - The LCP matrix
        b             - The LCP vector
        x0            - Initial guess
        max_iter      - Maximum number of iterations used to solve the LCP,
                        default to floor(problem_size / 2.0).
        tol_rel       - Relative tolerence, break if the error changes less 
                        than this between two iterations, otherwise it can
                        be configured to using gradient steps instead.
        tol_abs       - Absolute tolerence, break if the error drops below
                        this.
        subsolver     - Chooses which subsolver for solving the Newton
                        subsystem, these can be perturbation, zero, random
                        or approximation (not implemented).
        profile       - whether or not to profile while running, log 
                        convergence and perhaps other properties, number
                        of gradient steps, information about warm start etc.
        warmstart     - Class defining how to warm start the Fischer Newton
                        method. Should implement the entire warm start
                        procedure returning a new (and better) x0.
        gradient      - 
        searchmethod  - Defines the searchmethod used, the default Armijo
                        backtracking line search

    """

    def __init__(self):
        self.eps     = np.finfo(np.float64).eps
        self.rho     = np.finfo(np.float64).eps
        self.beta    = 0.001
        self.alpha   = 0.5

    def solve(self):
        pass

    def solve(self, A, b, x0, max_iter=None, tol_rel=0.0001, tol_abs=10*np.finfo(np.float64).eps, subsolver=None, profile=False, warmstart=None, gradient=None, searchmethod=None):

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

        if not A.shape[0] == A.shape[1]:
            raise NonSquareException("A is not a square matrix it has shape: {}".format(repr(A.shape)))

        if not A.shape[0] == x0.shape[0]:
            raise InvalidShapeException("A does not have the same dimensions as x0. A.shape: {}, while x0.shape: {}".format(repr(A.shape), repr(x0.shape)))

        if not b.shape[0] == x0.shape[0]:
            raise InvalidShapeException("x0 does not have the same dimensions as b. x0.shape: {}, while b.shape: {}".format(repr(x0.shape), repr(b.shape)))

        N = x0.shape[0]

        if max_iter:
            self.max_iter = max_iter
        else:
            self.max_iter = np.floor(N/2.0)

        if subsolver:
            self.subsolver = subsolver
        else:
            self.subsolver = SubsolverPerturbation()

        if gradient:
            self.use_gradient_steps = True
            self.take_gradient_step = False
        else:
            self.use_gradient_steps = False


        if searchmethod:
            self.searchmethod = searchmethod
        else:
            self.searchmethod = ArmijoLineSearch(merit_func=self.fischer)

        if warmstart:
            self.warm_start = True
            x0 = warmstart(A,b,x0, max_warm_iter)

        if profile:
            self.stat = StatStructure(max_iter+1, self.use_gradient_steps, self.warm_start)

        x = x0
        err = 1e20
        self.iterate = 1

        convergence = np.zeros(self.max_iter+1)

        while self.iterate <= self.max_iter:

            y = np.dot(A,x) + b

            phi = self.fischer(x,y)
            old_err = err
            err = 0.5*np.dot(phi.T,phi)

            if profile:
                convergence[self.iterate] = err

            if self.use_gradient_steps:
                self.take_gradient_step = False

            ##### Test the stopping criterias used #####
            rel_err = np.abs(err-old_err)/np.abs(old_err)

            if rel_err < tol_rel:
                if self.use_gradient_steps:
                    self.take_gradient_step = True
                else:
                    flag = 3
                    break

            if err < tol_abs:
                flag = 4
                break

            # Call the chosen subsolver
            (J, dx) = subsolver.solve(A, x, y, phi)

            nabla_phi = np.dot(phi.T, J)
            # Make it a column vector
            nabla_phi = nabla_phi.T

            # Perform gradient descent step if take_grad_step is defined
            if self.use_gradient_steps and self.take_gradient_step:
                self.searchmethod.alpha = gradient.grad_alpha
                dx = -nabla_phi
                if record_stat:
                    stat.gradient_steps +=1

            # Tests whether the search direction is below machine
            # precision.
            if np.max(np.abs(dx)) < self.eps:
                # If enabled and the relative change in error is low, use
                # a gradient descent step
                if take_gradient_step:
                    dx = -nabla_phi
                    grad_steps += 1
                else:
                    flag = 5
                    break

            # Test whether we are stuck in a local minima
            if np.linalg.norm(nabla_phi) < tol_abs:
                flag = 6
                break

            # Test whether our direction is a sufficient descent direction
            if np.dot(nabla_phi.T,dx) > -self.rho*(np.dot(dx.T, dx)):
                # Otherwise we should try gradient direction instead.
                if self.use_grad_steps:
                    dx = nabla_phi
                    stat.gradient_steps += 1
                    print ">>> Gradient descent step taken. (non descent)"
                else:
                    flag = 7
                    break

            # For testing sufficient decrease in function
            f_0 = err
            # Gradient of f
            grad_f = self.beta*np.dot(nabla_phi.T,dx)

            # Search method returns the new x^{(k+1)}, instead of
            # returning tau^{(k)}, which is more flexible?
            x[:] = self.searchmethod.findNextx(A, x, b, dx, f_0, grad_f)
            
            # Reset search alpha to original
            self.searchmethod.alpha = self.alpha

            self.iterate += 1

        if self.iterate >= self.max_iter:
            self.iterate -= 1
            flag = 8

        # Might want to change this to returning a
        # FischerNewtonSolution containing all the relevant stuff. Or
        # might this be too much of an overhead.
        return (x, err, self.iterate, flag, convergence[:self.iterate], msg[flag])


    def fischer(self,a,b):
        # Fischer function
        # Cite reference to it.
        return np.sqrt(np.power(a,2) + np.power(b,2)) - (a + b)
        # return np.sqrt(y**2+x**2) - y - x

    # def fischerLinearOperator(A,x,dx,phi,h,idx):
    #     fischer(np.dot(A[idx[0],:][:,idx[0]],(x[idx]+dx[idx]*h)),
    #             x[idx]+dx*h)-phi[idx]
        
    class SubsolverPerturbation:

        def __init__(self, gamma = 1e-28, eps=np.finfo(np.float64).eps):
            self.gamma = gamma
            self.eps   = eps

        def solve(self, A, x, y, phi):
            """Solves the sub Newton system of the Fischer Newton method
            in order to obtain the Newton direction dx."""

            # Bitmask of True/False values where the conditions are
            # satisfied, that is where phi and x are near
            # singular/singular.
            S = np.logical_and(np.abs(phi) < self.gamma,
                           np.abs(x) < self.gamma)

            N = np.size(x)
            
            px = np.zeros(x.shape)
            px[:] = x # Copy x
            assert px.shape == (N,1), 'px is not a column vector, it has shape: ' + \
                repr(px.shape)
            assert np.all( np.isreal(px) ), 'px is not real'

            direction = np.zeros(x.shape)
            direction[:] = np.sign(x)
            direction[ direction == 0 ] = 1
            px[S] = self.gamma * direction[S]

            p = np.zeros(px.shape)
            p = (px / (np.sqrt(np.power(y,2) + np.power(px,2))))-1
            assert p.shape == (N,1), 'p is not a column vector, it has shape: ' + \
                repr(p.shape)
            assert np.all(np.isreal(p)), 'p is not real'

            q = np.zeros(y.shape)
            q = (y / (np.sqrt(np.power(y,2) + np.power(px,2)))) - 1
            assert q.shape == (N,1), 'q is not a column vector, it has shape: ' + \
                repr(q.shape)
            assert np.all(np.isreal(q)), 'q is not real'
            
            J = np.dot(np.diag(p[:,0]),np.identity(N)) + np.dot(np.diag(q[:,0]),A)
            assert J.shape == (N,N), 'J is not a square matrix, it has shape: ' + \
                repr(J.shape)
            assert np.all(np.isreal(J)), 'J is not real'
            
            # Reciprocal conditioning number of J
            rcond = 1/np.linalg.cond(J)
            if np.isnan(rcond):
                rcond = np.nan_to_num(rcond)

            # Check conditioning of J, if bad use pinv to solve the system.
            if rcond < 2*self.eps:
                # dx = np.dot(np.linalg.pinv(J), (-phi))
                if np.any(np.isnan(J)):
                    print "J contains nan"
                if np.any(np.isnan(phi)):
                    print "phi contains nan"
                try:
                    dx = np.dot(spl.pinv(J), (-phi))
                except Exception:
                    print repr(dir(Exception))
                    dx = np.dot(spl.pinv2(J), (-phi))
            else:
                # if J.issparse():
                J = sps.csc_matrix(J)
                dx = spsl.dsolve.spsolve(J, (-phi), use_umfpack=True)
                J = J.toarray()
                # else:
                #     dx = np.linalg.solve(J,(-phi))

            dx = dx.reshape(dx.size,1)

            assert dx.shape == (N,1), 'dx is not a column vector, it has shape: ' + \
                repr(dx.shape)
            assert np.all(np.isreal(dx)), 'dx is not real'

            return (J, dx)
        
    class SubsolverRandom:
        pass

    class SubsolverZero:
        pass

    class SubsolverApproximation:
        # This has not yet been implemented.
        pass
