

class StatStructure():
    
    def __init__(self, size, gradient=False, warmstart=False):
        self.convergence = np.zeros()
        self.gradient_steps = 0
        self.warmstart_iter = 0
