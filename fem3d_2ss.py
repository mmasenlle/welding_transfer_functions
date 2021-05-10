import numpy as np
import control

# space state system from fem matrices (fenics3d_mesh.py)
class Fem3d_fenics:
    def __init__(self, path_data='data/'):
        try:
            self.X = np.load(path_data + 'XX3d.npy')
            self.A = np.load(path_data + 'AA3d.npy')
            self.B = np.load(path_data + 'BB3d.npy')
            self.Bv = np.load(path_data + 'BB3dv.npy')
            self.Teq = np.load(path_data + 'T03d.npy')
        except Exception as error:
            print('Fem3d_fenics::init->error: ', error)
    def get_C(self, xo, npts=4):
        C = np.zeros(self.X.shape[0])
        X_xo = self.X - np.array(xo)
        idx = np.argwhere(np.all((X_xo) == 0, axis=1))
        if idx.size:
            C[idx[0][0]] = 1
        else:
            dist = np.sqrt(X_xo[:,0]**2 + X_xo[:,1]**2 + X_xo[:,2]**2)
            # idx = np.argpartition(dist, 4)[:4]
            idx = np.argpartition(dist, npts)[:npts]
            distp = np.prod(dist[idx])
            C[idx] = (distp/np.sum(distp/dist[idx])) / dist[idx]
        return C
    def get_ss(self, xo):
        return control.ss(self.A, self.B, self.get_C(xo), 0)
    def get_ss_v(self, xo):
        return control.ss(self.A, self.Bv, self.get_C(xo), 0)
