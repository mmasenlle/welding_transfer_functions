import numpy as np

class Filter1:
    def __init__(self, b, a):
        self.b = np.array(b)
        self.a = np.array(a)
        self.x = np.zeros(self.b.size - 1)
    def step(self, u):
        y = (self.b[0] * u) + self.x[0]
        for i in range(1, self.x.size):
            self.x[i-1] = (self.b[i] * u) + self.x[i] - (self.a[i] * y)
        self.x[-1] = (self.b[-1] * u) - (self.a[-1] * y)
        return y