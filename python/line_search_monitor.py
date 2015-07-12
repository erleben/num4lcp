import numpy as np

class LineSearchMonitor(object):
    def __init__( self, alpha=.5, beta=.0001, tau_min=np.finfo(float).eps ):
        self.alpha            = alpha
        self.beta             = beta
        self.tau_min          = tau_min
        self.line_search_iter = 0
        self.tau              = 1.0
        self.fk               = np.finfo(np.float64).max
        self.f0               = np.finfo(np.float64).max
        self.grad_f0          = np.finfo(np.float64).max

    def init( self, f0, grad_f0 ):
        self.f0 = f0
        self.grad_f0 = self.beta*grad_f0

    def reset( self ):
        self.line_search_steps = 0
        self.tau = 1.0
        
    def finished( self, fk ):
        self.fk = fk
        return self.step_length_found() or self.check_step_length()

    def step_length_found( self ):
        if self.fk <= self.f0 + self.tau*self.grad_f0:
            return True
        else:
            self.tau*=self.alpha
            return False

    def check_step_length( self ):
        if self.tau*self.tau < self.tau_min:
            return True
        else:
            return False

    def step( self ):
        self.line_search_iter += 1

    def step_length( self ):
        return self.tau

    def line_search_steps( self ):
        return self.line_search_iter
