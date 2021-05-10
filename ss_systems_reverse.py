
import fem3d_2ss
import numpy as np
fem2ss = fem3d_2ss.Fem3d_fenics('test_open/')

p1=np.load('test_open/p1.npy')
p2=np.load('test_open/p2.npy')
# p1r=np.load('test_open/p1r.npy')
vars={}
vars['A'] = fem2ss.A
vars['B'] = fem2ss.B
vars['Bv'] = fem2ss.Bv
vars['C1'] = fem2ss.get_C(p1)
vars['C2'] = fem2ss.get_C(p2, 2)
# vars['C1r'] = fem2ss.get_C(p1r)
# vars['Ar'] = np.load('test_open/AA3dr.npy')
# vars['Br'] = np.load('test_open/BB3dr.npy')
# vars['Bvr'] = np.load('test_open/BB3dvr.npy')
vars['od'] = np.load('test_open/open_loop_data.npy')
vars['T0'] = np.load('test_open/T03d.npy')
vars['T0r'] = np.load('test_open/T03dr.npy')

# X_xo = fem2ss.X - np.array(p2)
# idx = np.argwhere(np.all((X_xo) == 0, axis=1))
#
# dist = np.sqrt(X_xo[:,0]**2 + X_xo[:,1]**2 + X_xo[:,2]**2)
# idx = np.argpartition(dist, 2)[:2]
# fem2ss.X[idx]
# fem2ss.Teq[idx[:]]
#
# distp = np.prod(dist[idx])
# C = np.zeros(fem2ss.X.shape[0])
# C[idx] = (distp/np.sum(distp/dist[idx])) / dist[idx]
# C @ fem2ss.Teq
# C @ vars['T0']
# C @ vars['T0r']
# vars['C1'] @ vars['T0']
# vars['C1'] @ vars['T0r']
# vars['C2'] @ vars['T0']
# vars['C2'] @ vars['T0r']
# vars['od'][-2:]


import scipy.io
scipy.io.savemat('ss_vars_reverse.mat', vars)
