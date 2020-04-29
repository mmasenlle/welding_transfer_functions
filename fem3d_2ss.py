import numpy as np
import control

# space state system from fem matrices (fenics3d_mesh.py)
class Fem3d_fenics:
    def __init__(self):
        self.X = np.load('XX3d.npy')
        self.A = np.load('AA3d.npy')
        self.B = np.load('BB3d.npy')
        self.Bv = np.load('BB3dv.npy')
    def get_ss(self, xo):
        C = np.zeros(self.X.shape[0])
        C[np.argwhere(np.all((self.X-np.array(xo))==0, axis=1))[0][0]] = 1
        return control.ss(self.A, self.B, C, 0)
    def get_ss_v(self, xo):
        C = np.zeros(self.X.shape[0])
        C[np.argwhere(np.all((self.X-np.array(xo))==0, axis=1))[0][0]] = 1
        return control.ss(self.A, self.Bv, C, 0)
