import numpy as np
import control

# space state system from fem matrices (fenics2d_ss.py)
class Fem2d_fenics:
    def __init__(self):
        self.X = np.load('data/2d/XX2d.npy')
        self.A = np.load('data/2d/AA2d.npy')
        self.B = np.load('data/2d/BB2d.npy')
        self.Bv = np.load('data/2d/BB2dv.npy')
    def get_C(self, xo):
        C = np.zeros(self.X.shape[0])
        X_xo = self.X - np.array(xo)
        idx = np.argwhere(np.all((X_xo) == 0, axis=1))
        if idx.size:
            C[idx[0][0]] = 1
        else:
            dist = np.sqrt(X_xo[:,0]**2 + X_xo[:,1]**2)
            idx = np.argpartition(dist, 4)[:4]
            distp = np.prod(dist[idx])
            C[idx] = (distp/np.sum(distp/dist[idx])) / dist[idx]
        return C
    def get_ss(self, xo):
        return control.ss(self.A, self.B, self.get_C(xo), 0)
    def get_ss_v(self, xo):
        return control.ss(self.A, self.Bv, self.get_C(xo), 0)
