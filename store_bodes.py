import numpy as np
import control


import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/16x8x2/')
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/')

inp = np.array((.02, .02, 0))
pts = ((0,0.003,0),(0,0.005,0),(0,0.010,0),(0.005,0,0),(0.010,0,0),(0.020,0,0))
w = np.logspace(-3,1,100)

legends=[]
bodes_data=[]
for p in pts[:2]:
    legends.append(str(p))
    bodes_data.append(np.array(control.bode(fem2ss.get_ss(inp + p), w, dB=True, Plot=False)))
np.save('bodes_data',bodes_data)
exit(0)
# show plots
bds = np.load('bodes_data.npy')
import matplotlib.pyplot as plt
fig = plt.figure('Bode T/Pot')
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.set_title('Bode T/Pot')
for bd in bds:
    ax1.semilogx(bd[2], 20*np.log10(bd[0]))
    ax2.semilogx(bd[2], bd[1]*(180/np.pi))
ax2.legend(legends)
ax1.set_xlim((w[0],w[-1]))
ax2.set_xlim((w[0],w[-1]))
