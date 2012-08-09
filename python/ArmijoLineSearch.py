import numpy as np

"""
Copyright 2012, Michael Andersen, DIKU. michael (at) diku (dot) dk
"""

class ArmijoLineSearch():

    # def __init__():
    #     pass

    def __init__(self, alpha=0.5, gamma=1e-28, tau=1.0, merit_func=None):
        self.alpha = alpha
        self.gamma = gamma
        self.tau   = tau
        self.merit_func = merit_func

    def findNextx(self, A, x, b, dx, f_0, grad_f):
        
        tau_k = self.tau

        x_k = x[:]

        while True:
            x_k = np.maximum(0, x + dx*tau_k)
            y_k = np.dot(A,x_k)+b

            # Calculate phi value
            phi_k = self.merit_func(x_k,y_k)

            # Calculate the natural merit function
            f_k = 0.5*(np.dot(phi_k.T,phi_k))

            # Test Armijo condition for sufficient decrease
            if np.all(f_k <= f_0 + tau_k*grad_f):
                break
            
            # Test whether the stepsize has become too small
            if tau_k*tau_k < self.gamma:
                break

            tau_k *= self.alpha
            
        return x_k
