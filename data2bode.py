import numpy as np
import matplotlib.pyplot as plt
import control


k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)

ops=np.load('data/output_points0.npy')

X=np.load('data/output_data_pow.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,int(np.log10(N)),100).astype(int)
w = wf[idx]
s = 1j*w
ft_pot = lambda x,y,z: 1/(rho*cp)*np.exp(-(v*(-x)+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))

import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics('data/v2/16x8x2/')
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v2/16x8x2/')

U = np.fft.fft(X[:,1])
for i in range(3,X.shape[1]):
    H1 = U * np.fft.fft(X[:,i]) / U**2
    plt.figure('Bode T/Pow ' + str(ops[i - 3]))
    control.bode_plot((control.frd(H1[idx], w), control.frd(ft_pot(*ops[i - 3]), w), fem2ss.get_ss(ops[i - 3]+(.02,.02,0))), w, dB=True)
    plt.legend(('fem', 'fdt', 'ss'))


U = np.fft.fft(X[:,1])
H1 = U * np.fft.fft(X[:,7]) / U**2
H2 = U * np.fft.fft(X[:,8]) / U**2
plt.figure('Bode T/Pow simetricity comparative 1')
control.bode_plot((control.frd(H1[idx], w), control.frd(ft_pot(*ops[4]), w), #fem2ss.get_ss(ops[4]+(.02,.02,0)),
                   control.frd(H2[idx], w), control.frd(ft_pot(*ops[5]), w)), #fem2ss.get_ss(ops[5]+(.02,.02,0))),
                  w, dB=True)
# plt.legend(('fem +5', 'fdt +5', 'ss +5', 'fem -5', 'fdt -5', 'ss -5'))
plt.legend(('fem +5', 'fdt +5', 'fem -5', 'fdt -5'))

n=4
ops_all=np.load('output_points_0.npy')-(.02,.02,0)
X_all=np.load('output_data_pow_0.npy')
for i in range(1,n):
    ops_all = np.vstack((ops_all, np.load('output_points_'+str(i)+'.npy') - (.02, .02, 0)))
    X_all = np.hstack((X_all, np.load('output_data_pow_'+str(i)+'.npy')[:,3:]))



U = np.fft.fft(X_all[:,1])
H1 = U * np.fft.fft(X_all[:,6]) / U**2
H2 = U * np.fft.fft(X_all[:,8]) / U**2
plt.figure('Bode T/Pow simetricity comparative 2')
control.bode_plot((control.frd(H1[idx], w), control.frd(ft_pot(*ops_all[3]), w), #fem2ss.get_ss(ops[4]+(.02,.02,0)),
                   control.frd(H2[idx], w), control.frd(ft_pot(*ops_all[5]), w)), #fem2ss.get_ss(ops[5]+(.02,.02,0))),
                  w, dB=True)
# plt.legend(('fem +5', 'fdt +5', 'ss +5', 'fem -5', 'fdt -5', 'ss -5'))
plt.legend(('fem +5', 'fdt +5', 'fem -5', 'fdt -5'))
# plt.xlim((0.01,10))
# plt.ylim((-50,0))
# plt.axis('equal')

