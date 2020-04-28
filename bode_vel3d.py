import numpy as np
import matplotlib.pyplot as plt
import control

k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)
q=100

'''
ops=np.load('output_points0.npy')

X=np.load('output_data_vel.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,3,100).astype(int)
w = wf[idx]
'''
w = np.logspace(-3,2,100)
s = 1j*w


ft_vel = lambda x,y,z: -q/(rho*cp)*(v*np.sqrt(x**2 + y**2 + z**2) + (-x)*np.sqrt(4*a*s + v**2))*np.exp(-(v*(-x) + np.sqrt(4*a*s + v**2)*np.sqrt(x**2 + y**2 + z**2))/(2*a))/(4*np.pi*a**2*np.sqrt(4*a*s + v**2)*np.sqrt(x**2 + y**2 + z**2))

R = lambda x,y,z: np.sqrt(x**2+y**2+z**2)
D = lambda s: np.sqrt(4*a*s + v**2)
XX = lambda x,y,z,s: (v*R(x,y,z) - x*D(s))*np.exp(-R(x,y,z)*D(s)/(2*a))/D(s)
ft_vel2 = lambda x,y,z: -q/(rho*cp)*np.exp(v*x/(2*a))/(4*s*np.pi*a**2*R(x,y,z))*(XX(x,y,z,s)-XX(x,y,z,0))

ops = ((-.01,0,0),(.01,0,0),(.02,0,0),(.04,0,0),(0,-.01,0),(0,.01,0))

#U = np.fft.fft(X[:,2])
for i in range(len(ops)):
    #H1 = U * np.fft.fft(X[:,i]) / U**2
    plt.figure('Bode T/Vel ' + str(ops[i]))
    control.bode_plot((control.frd(ft_vel2(*ops[i]), w),control.frd(ft_vel(*ops[i]), w)),w, dB=True)
    plt.legend(('fdt2', 'fdt'))
    plt.show()
