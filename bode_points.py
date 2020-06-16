import numpy as np
import matplotlib.pyplot as plt
import control


import fem3d_2ss
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/16x8x2/')
# fem2ss_32 = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/')

inp = np.array((.02, .02, 0))
pts = ((0,0.003,0),(0,0.005,0),(0,0.010,0),(0.005,0,0),(0.010,0,0),(0.020,0,0))
w = np.logspace(-3,1,100)

k=24.0
rho=7925
cp=460
v=0.005
a=k/(rho*cp)
s = 1j*w
ft_pot = lambda x,y,z: 1/(rho*cp)*np.exp(-(v*(-x)+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))

pts = ((0,0.003,0),(0,-0.003,0),(0,0.005,0),(0,-0.005,0),(0,0.010,0),(0,-0.010,0))
fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/irregular/')
fem2ss_32 = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/1500/')
for p in pts:
    plt.figure('Bode Comp ' + str(p))
    control.bode_plot((fem2ss.get_ss(inp + p), fem2ss_32.get_ss(inp + p), control.frd(ft_pot(*p), w)), w, dB=True)
    plt.legend(('irregular', '32', 'ft'))


# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/1500/')

fem2ss = fem3d_2ss.Fem3d_fenics('./')
pts = ((0,0.003,0),(0,-0.003,0),(0,0.005,0),(0,-0.005,0),(0,0.010,0),(0,-0.010,0))

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
vars['C_00_03'] = fem2ss.get_C((.02,.023,0))
vars['C_00__03'] = fem2ss.get_C((.02,.017,0))
vars['C_00_05'] = fem2ss.get_C((.02,.025,0))
vars['C_00__05'] = fem2ss.get_C((.02,.015,0))
vars['C_00_07'] = fem2ss.get_C((.02,.027,0))
vars['C_00__07'] = fem2ss.get_C((.02,.013,0))
vars['C_00_10'] = fem2ss.get_C((.02,.03,0))
vars['C_00__10'] = fem2ss.get_C((.02,.01,0))
vars['C_00_13'] = fem2ss.get_C((.02,.033,0))
vars['C_00__13'] = fem2ss.get_C((.02,.007,0))
vars['C_03_00'] = fem2ss.get_C((.023,.02,0))
vars['C__03_00'] = fem2ss.get_C((.017,.02,0))
vars['C_05_00'] = fem2ss.get_C((.025,.02,0))
vars['C__05_00'] = fem2ss.get_C((.015,.02,0))
vars['C_07_00'] = fem2ss.get_C((.027,.02,0))
vars['C__07_00'] = fem2ss.get_C((.013,.02,0))
vars['C_10_00'] = fem2ss.get_C((.03,.02,0))
vars['C__10_00'] = fem2ss.get_C((.01,.02,0))
vars['C_12_00'] = fem2ss.get_C((.032,.02,0))
vars['C__12_00'] = fem2ss.get_C((.008,.02,0))
vars['C_15_00'] = fem2ss.get_C((.035,.02,0))
vars['C_17_00'] = fem2ss.get_C((.037,.02,0))
vars['C_20_00'] = fem2ss.get_C((.04,.02,0))
vars['C_25_00'] = fem2ss.get_C((.045,.02,0))
vars['C_30_00'] = fem2ss.get_C((.05,.02,0))
vars['C_40_00'] = fem2ss.get_C((.06,.02,0))
import scipy.io
scipy.io.savemat('ss_vars.mat', vars)

bodes_data=[]
for p in pts:
    bodes_data.append(control.bode(fem2ss.get_ss(inp + p), w, dB=True, Plot=False))