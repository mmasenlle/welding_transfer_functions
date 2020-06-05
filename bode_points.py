import numpy as np
import matplotlib.pyplot as plt
import control


import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/16x8x2/')
fem2ss_32 = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/')

inp = np.array((.02, .02, 0))
pts = ((0,0.003,0),(0,0.005,0),(0,0.010,0),(0.005,0,0),(0.010,0,0),(0.020,0,0))
w = np.logspace(-3,1,100)

for p in pts[:2]:
    plt.figure('Bode Comp ' + str(p))
    control.bode_plot((fem2ss.get_ss(inp + p), fem2ss_32.get_ss(inp + p)), w, dB=True)
    plt.legend(('16', '32'))


fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/1500/')

systemsQ = list()
systemsV = list()
legends = list()

for p in pts:
    systemsQ.append(fem2ss.get_ss(inp + p))
    systemsV.append(fem2ss.get_ss_v(inp + p))
    legends.append(str(p))


plt.figure('Bode T/Pot')
control.bode_plot(systemsQ, w, dB=True)
plt.legend(legends)
plt.figure('Bode T/Vel')
control.bode_plot(systemsV, w, dB=True)
plt.legend(legends)

# fem2ss=fem2ss_32
vars={}
vars['A'] = fem2ss.A
vars['B'] = fem2ss.B
vars['Bv'] = fem2ss.Bv
vars['C_00_05'] = fem2ss.get_C((.02,.025,0))
vars['C_00_07'] = fem2ss.get_C((.02,.027,0))
vars['C_00_10'] = fem2ss.get_C((.02,.03,0))
vars['C_05_00'] = fem2ss.get_C((.025,.02,0))
vars['C_10_00'] = fem2ss.get_C((.03,.02,0))
vars['C_20_00'] = fem2ss.get_C((.04,.02,0))
import scipy.io
scipy.io.savemat('ss_vars.mat', vars)

bodes_data=[]
for p in pts:
    bodes_data.append(control.bode(fem2ss.get_ss(inp + p), w, dB=True, Plot=False))