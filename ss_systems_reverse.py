
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

vars['AAA'] = np.load('test_open/AAA.npy')
vars['BB'] = np.load('test_open/BB.npy')
vars['BBv'] = np.load('test_open/BBv.npy')
vars['TT0'] = np.load('test_open/TT0.npy')
vars['TT0r'] = np.load('test_open/TT0r.npy')


import scipy.io
# scipy.io.savemat('ss_vars_reverse.mat', vars)
scipy.io.savemat('ss_plants2.mat', vars)


vars={}
import numpy as np
vars['cd'] = np.load('ctrl_qft_manu4.npy')
import scipy.io
scipy.io.savemat('ctrl_qft_manu5.mat', vars)

# vars={}
# vars['step10'] = np.load('test_open/ctrl_steps_data.npy')
# vars['step100'] = np.load('test_open/ctrl_steps_data100.npy')
# import scipy.io
# scipy.io.savemat('ctrl_steps.mat', vars)