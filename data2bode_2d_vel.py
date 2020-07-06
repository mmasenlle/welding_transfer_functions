import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv
import control

# T0=np.load('T0_2d.npy')
# X = np.load('X_2d.npy')
# nds1f = np.argwhere((X[:, 1] == 0.02))
# plt.plot(X[nds1f,0],T0[nds1f])
# nds1f = np.argwhere((X[:, 0] == 0.02))
# plt.plot(X[nds1f,1],T0[nds1f])


k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)
h=1000
b=2*h/(rho*cp)

ops=np.load('data/output_points_2d_vel.npy')

X=np.load('data/output_data_vel_2d.npy')

# plt.plot(X[:10000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
# idx = np.logspace(0,int(np.log10(N)),100).astype(int)
idx = np.logspace(0,4,100).astype(int)
idx[0] = 0
w = wf[idx]
s = 1j*w
s[0] = 1j*1e-10
K=k
q=100000

R = lambda x,y: np.sqrt(x**2+y**2)
D = lambda s: np.sqrt(4*a*s + 4*a*b + v**2)

chi = lambda x,y,s: v*kv(0,R(x,y)*D(s)/(2*a)) - (x/R(x,y))*(2*np.abs(v)-(v**2)/D(s))*kv(1,R(x,y)*D(s)/(2*a))
ft_vel15 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(chi(x,y,s) - chi(x,y,0))

import fem2d_2ss
fem2ss = fem2d_2ss.Fem2d_fenics()

U = np.fft.fft(X[:,2])
for i in range(3,X.shape[1]):
    H1 = U * np.fft.fft(X[:,i]) / U**2
    plt.figure('Bode T/V ' + str(ops[i - 3]))
    plt.title('Bode T/V output point=' + str(ops[i - 3]))
    control.bode_plot((control.frd(H1[idx], w), control.frd(ft_vel15(*ops[i - 3]), w), fem2ss.get_ss_v(ops[i - 3]+(.02,.02))),w, dB=True)
    plt.legend(('fem', 'fdt', 'ss'))

if False:
    vars={}
    vars['A'] = fem2ss.A
    vars['B'] = fem2ss.B
    vars['Bv'] = fem2ss.Bv
    vars['C_00_03'] = fem2ss.get_C((.02,.023))
    vars['C_00__03'] = fem2ss.get_C((.02,.017))
    vars['C_00_05'] = fem2ss.get_C((.02,.025))
    vars['C_00__05'] = fem2ss.get_C((.02,.015))
    vars['C_00_07'] = fem2ss.get_C((.02,.027))
    vars['C_00__07'] = fem2ss.get_C((.02,.013))
    vars['C_00_10'] = fem2ss.get_C((.02,.03))
    vars['C_00__10'] = fem2ss.get_C((.02,.01))
    vars['C_00_13'] = fem2ss.get_C((.02,.033))
    vars['C_00__13'] = fem2ss.get_C((.02,.007))
    vars['C_03_00'] = fem2ss.get_C((.023,.02))
    vars['C__03_00'] = fem2ss.get_C((.017,.02))
    vars['C_05_00'] = fem2ss.get_C((.025,.02))
    vars['C__05_00'] = fem2ss.get_C((.015,.02))
    vars['C_07_00'] = fem2ss.get_C((.027,.02))
    vars['C__07_00'] = fem2ss.get_C((.013,.02))
    vars['C_10_00'] = fem2ss.get_C((.03,.02))
    vars['C__10_00'] = fem2ss.get_C((.01,.02))
    vars['C_12_00'] = fem2ss.get_C((.032,.02))
    vars['C__12_00'] = fem2ss.get_C((.008,.02))
    vars['C_15_00'] = fem2ss.get_C((.035,.02))
    vars['C_17_00'] = fem2ss.get_C((.037,.02))
    vars['C_20_00'] = fem2ss.get_C((.04,.02))
    vars['C_25_00'] = fem2ss.get_C((.045,.02))
    vars['C_30_00'] = fem2ss.get_C((.05,.02))
    vars['C_40_00'] = fem2ss.get_C((.06,.02))
    import scipy.io
    scipy.io.savemat('ss_vars_2d.mat', vars)
