import numpy as np
from line_search_monitor import LineSearchMonitor

class NonmonotoneLineSearchMonitor(LineSearchMonitor):
    def __init__( self, alpha=.5, beta=.0001, memory_size=8, pre_newton_steps=5, tau_min=np.finfo(float).eps ):
        super().__init__( alpha, beta, tau_min )
        # self.memory = np.array(np.zeros(memory))
        self.memory_size = memory_size
        self.pre_newton_steps = pre_newton_steps
        self.outer_iter = 1 # This could be passed from outside on reset
        
        
    def init( self, f0, grad_f0 ):
        super().init(f0, grad_f0)
        if self.outer_iter == 1:
            self.init_memory();

    def init_memory( self ):
        self.memory = np.array([self.f0]*self.memory_size)

    def reset( self ):
        super().reset()
        self.outer_iter += 1
        self.memory[self.outer_iter % self.memory_size] = self.f0
        if self.outer_iter > pre_newton_steps:
            self.f_max = np.max(self.memory)
        else:
            self.f_max = self.f0
