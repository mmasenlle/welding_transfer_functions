
import control
import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics('data_isotermas/box2/')


sys1200 = ((fem2ss.get_ss((0.034018,0.02,0)),fem2ss.get_ss((0.02,0.022354,0))),
           (fem2ss.get_ss_v((0.034018,0.02,0)),fem2ss.get_ss_v((0.02,0.022354,0))))

control.bode(sys1200)



vars['C_1200y'] = fem2ss.get_C((0.02,0.022354,0))
vars['C_1000x'] = fem2ss.get_C((0.037608,0.02,0))
vars['C_1000y'] = fem2ss.get_C((0.02,0.022464,0))
vars['C_800x'] = fem2ss.get_C((0.041366,0.02,0))
vars['C_800y'] = fem2ss.get_C((0.02,0.022574,0))
vars['C_600x'] = fem2ss.get_C((0.048639,0.02,0))
vars['C_600y'] = fem2ss.get_C((0.02,0.022834,0))


import scipy.io
scipy.io.savemat('ss_vars_isotherms.mat', vars)
