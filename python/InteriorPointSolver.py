import numpy as np
from StatStructure import StatStructure

"""
Copyright 2012, Michael Andersen, DIKU. michael (at) diku (dot) dk
"""


class InteriorPointSolver:

    def __init__(self):
        self.eps     = np.finfo(np.float64).eps
    # Human readable dictionary of exit messages
        self.msg = {1 : 'preprocessing',  # flag = 1
                    2 : 'iterating',      # flag = 2
                    3 : 'relative',       # flag = 3
                    4 : 'absolute',       # flag = 4
                    5 : 'stagnation',     # flag = 5
                    6 : 'local minima',   # flag = 6
                    7 : 'nondescent',     # flag = 7
                    8 : 'maxlimit',       # flag = 8
                    }


    def solve(self):
        pass

    def solve(self, A, b, x0, max_iter=0, tol_rel=0.00001, tol_abs=np.finfo(np.float64).eps*10, subsolver=None, profile=False, warmstart=None, gradient=None, searchmethod=None ):

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
            # max_iter = problem_size seems to work fine
            self.max_iter = N

        if subsolver:
            self.subsolver = subsolver
        else:
            self.subsolver = self.CentralTrajectory()

        if gradient:
            self.use_gradient_steps = gradient
            self.take_gradient_step = False
        else:
            self.use_gradient_steps = gradient

        # if searchmethod:
        #     self.searchmethod = searchmethod
        # else:
        #     # This is not correct!!!
        #     self.searchmethod = ArmijoLineSearch(merit_func=self.fischer)

        if warmstart:
            self.warm_start = True
            x0 = warmstart(A,b,x0, max_warm_iter)
        else:
            self.warm_start = False

        if profile:
            self.stat = StatStructure(max_iter+1, self.use_gradient_steps, self.warm_start)
        
        # Setup the things we need to compute the solution.
        y0 = np.dot(A,x0)+b
        e  = np.ones((N,1))
        I  = np.identity(N)

        # Make sure initial values are feasible that is x0, y0 \in S_{++}?
        if np.dot(x0.T,y0) <= 0:
            W = np.diag(y0[:,0])
            M = np.diag(x0[:,0])
            # Jtmp1 = [A | -I]
            Jtmp1 = np.hstack((A,-I))
            # Jtmp2 = [W | M]
            Jtmp2 = np.hstack((W,M))
            #     | Jtmp1 |
            # J = | ----- |
            #     | Jtmp2 |
            J = np.vstack((Jtmp1,Jtmp2))
            #      | A*x0-y0+b |
            # Fk = | --------- |
            #      | M*W*e - 1 |
            Fk = np.vstack((np.dot(A,x0) - y0 + b,np.dot(np.dot(M,W),e)-1))

            # Inverse condition number of J
            rcond = 1/np.linalg.cond(J)

            # Is the problem well conditioned, if not use pinv.
            if rcond < 2*self.eps:
                p = np.dot(-np.linalg.pinv(J),Fk)
            else:
                p = -np.linalg.solve(J,Fk)

            dx = p[:N]
            dy = p[N:]

            x0 = np.maximum(1, np.abs(x0 + dx))
            y0 = np.maximum(1, np.abs(y0 + dy))

        old_err = 1e20

        # Copy over x0 and y0
        x = x0[:]
        y = y0[:]

        flag = 2

        self.iterate = 1

        while self.iterate < self.max_iter: # We start with 1 not zero as we
                                   # are supposed to.

            # Merit value, complementarity measure, error.
            err = np.dot(x.T, y)

            if err <= 0:
                print ">>> ERROR x'*y = %2.5f" % (err)
                return

            if profile:
                self.stat.convergence[self.iterate - 1] = err

            if np.abs(err - old_err)/np.abs(old_err) < tol_rel:
                if profile:
                    self.stat.flag = 3
                break

            if err < tol_abs:
                if profile:
                    self.stat.flag = 4
                break
        
            sigma = np.float64(1)/self.iterate # Goes towards zero as we iterate
            gamma = sigma*np.dot(x.T,y)/N # Should go towards zero as we iterate

            ##### Setup and solve the Newton system ##########

            p = self.subsolver.solve(A,b,x,y, sigma, gamma)
        
            # ##### Setup the Newton system ###############

            # # if 'central_trajectory' == str.lower(solver):
            # W = np.diag(y[:,0])
            # M = np.diag(x[:,0])

            # Jtmp1 = np.hstack((A,-I))
            # Jtmp2 = np.hstack((W, M))
            # J = np.vstack((Jtmp1,Jtmp2))
            # Fk = np.vstack((np.dot(A,x) - y + b,
            #                 np.dot(np.dot(M,W),e)-gamma*e))
        
            # ##### Solve the Newton system ###############

            # p = -np.linalg.solve(J,Fk)
            # # linalg.inv is very slow.
            # # p = np.dot(-np.linalg.inv(J), Fk)

            # # elif 'schur_complement':

            # # Needs to be implmented

            # # elif 'potential_reduction':

            # # Needs to be implmented

            # # else:
            # #     print ">>> ERROR: Unknown subsolver: {}".format(solver)

            dx = p[:N]
            dy = p[N:]

            # Controls how much we step back from the maximum step
            # length. Such that we still respect the conditions
            # x^{(k+1)}_i > and y^{(k+1)}_i > 0. 0 < alpha < 1
            alpha = np.float64(self.iterate)/(self.max_iter+1.0) # Goes towards one as we iterate

            xx = ((1-alpha)*x - x)/dx
            yy = ((1-alpha)*y - y)/dy

            if not np.size(xx[dx < 0]) == 0:
                taux = np.minimum( 1, np.min(xx[dx < 0]))
            else:
                taux = 1

            if not np.size(yy[dy < 0]) == 0:
                tauy = np.minimum(1, np.min(yy[dy < 0]))
            else:
                tauy = 1

            # Verify that the step length is valid.
            tau = np.min([taux, tauy])

            if tau < 0 or tau > 1:
                print ">>> ERROR tau = %2.20f" % tau
                return
        
            ### Take the Newton step ###############
            x = x + taux*dx
            y = y + tauy*dy

            # Verify that the x and y vector are still valid. x > 0 and y
            # > 0.
            if np.min(x) <= 0:
                print ">>> ERROR Min x = %2.20f" % (np.min(x))
                return

            if np.min(y) <= 0:
                print ">>> ERROR Min y = %2.20f" % (np.min(y))
                return

            self.iterate += 1
    
        if self.iterate >= self.max_iter:
            self.iterate -= 1
            if profile:
                self.stat.flag = 8
        
        if profile:
            self.stat.msg = self.msg[self.stat.flag]
            self.stat.iterate = self.iterate

        # Return the solution, err, iterates, exit flag, convergence and
        # human readable message.
        if profile:
            return (x, err, self.stat)
        else:
            return (x, err)

    class CentralTrajectory:
        
        def __init__(self):
            pass

        def solve(self):
            pass

        def solve(self, A, b, x, y, sigma, gamma):
            # We initialize these each time we call the method, might
            # be very stupid.
            N = A.shape[0]
            e = np.ones((N,1))

            W = np.diag(y[:,0])
            M = np.diag(x[:,0])

            Jtmp1 = np.hstack((A,-np.identity(N)))
            Jtmp2 = np.hstack((W, M))
            J = np.vstack((Jtmp1,Jtmp2))
            Fk = np.vstack((np.dot(A,x) - y + b,
                            np.dot(np.dot(M,W),e)-gamma*e))

            # Should we use rcond?

            p = -np.linalg.solve(J,Fk)
        
            return p
