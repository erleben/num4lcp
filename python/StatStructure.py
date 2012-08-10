import numpy as np

class StatStructure():
    
    def __init__(self, size, gradient=False, warmstart=False):
        self.convergence = np.zeros(size)
        self.flag = 0
        self.msg = ''
        self.iterate = 0

        if gradient:
            self.gradient = gradient
            self.gradient_steps = {'non descent': 0,
                                   'search direction': 0,
                                   'gradient descent': 0}
        else:
            self.gradient = gradient
        if warmstart:
            self.warmstart = warmstart
            self.warmstart_iter = 0
        else:
            self.warmstart = warmstart

    def __repr__(self):
        gradient_string = ''
        warmstart_string = ''
        if self.gradient:
            gradient_string = 'gradient steps: '
            for key in self.gradient_steps:
                if self.gradient_steps[key] != 0:
                    gradient_string += key + repr(self.gradient_steps[key]) + ' '
        if self.warmstart:
            warmstart_string = ', warm start iterations: {}'.format(self.warmstart_iter)
        return 'Iteration: {}, message: {}, {} {}'.format(self.iterate, self.msg, gradient_string, warmstart_string)

