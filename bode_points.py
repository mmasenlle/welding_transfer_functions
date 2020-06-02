import numpy as np
import matplotlib.pyplot as plt
import control


import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics()

inp = np.array((.02, .02, 0))
pts = ((0,0.003,0),(0,0.005,0),(0,0.007,0),(0,0.010,0),(0.003,0,0),(0.005,0,0),(0.007,0,0),(0.010,0,0))
systemsQ = list()
systemsV = list()
legends = list()

for p in pts:
    systemsQ.append(fem2ss.get_ss(inp + p))
    systemsV.append(fem2ss.get_ss_v(inp + p))
    legends.append(str(p))

w = np.logspace(-3,2,100)
plt.figure('Bode T/Pot')
control.bode_plot(systemsQ, w, dB=True)
plt.legend(legends)
plt.figure('Bode T/Vel')
control.bode_plot(systemsV, w, dB=True)
plt.legend(legends)

