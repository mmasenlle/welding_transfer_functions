import numpy as np
import matplotlib.pyplot as plt
import control

k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)
q=1000
dl=0.001

ops=np.load('output_points0.npy')

X=np.load('output_data_pow.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,3,100).astype(int)
w = wf[idx]
s = 1j*w
ft_pot = lambda x,y,z: 1/(rho*cp)*np.exp(-(v*x+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))

Mf=None
af=None
frds_fdt=[]
frds_fem=[]

U = np.fft.fft(X[:,1])
for i in range(3,X.shape[1]):
    H1 = U * np.fft.fft(X[:,i]) / U**2
    plt.figure('Bode T/Pow ' + str(ops[i - 3]))
    control.bode_plot((control.frd(H1[idx], w),control.frd(ft_pot(*ops[i - 3]), w)),w, dB=True)
    plt.legend(('fem', 'fdt'))
    # frds_fem.append(control.frd(H1[:N], wf))
    # frds_fdt.append(control.frd(ft_pot(*ops[i - 3]), wf))
    #
    # H1 = H1[:N]
    #
    # Mf = np.vstack((Mf, 20*np.log10(abs(H1)))) if Mf is not None else 20*np.log10(abs(H1))
    # # Mf = np.vstack((Mf, abs(H1))) if Mf is not None else abs(H1)
    # af = np.vstack((af, np.angle(H1)*180/np.pi)) if af is not None else np.angle(H1)*180/np.pi


plt.figure('Magnitude T/Pow')
for i in range(X.shape[1]-3):
    plt.semilogx(wf, Mf[i,:])
plt.legend(('x=-.01','x=.01','x=.02','x=.04','y=-.005','y=.005','y=.01','z=.005'))

k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)
q=1000
dl=0.001

w = np.logspace(-3,1,100)
s = 1j*w

# ft_pot = lambda x,y,z: q/(rho*cp)*np.exp(-(v*x+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))
ft_pot = lambda x,y,z: 1/(rho*cp)*np.exp(-(v*x+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))

ops=np.load('output_points0.npy')[:2]
for i,op in enumerate(ops):
    plt.figure('Magnitude T/Pow '+str(op))
    fr = ft_pot(*op)
    plt.subplot(211)
    plt.semilogx(wf, Mf[i, :])
    plt.semilogx(w, 20 * np.log10(abs(fr)))
    plt.subplot(212)
    plt.semilogx(wf, af[i, :])
    plt.semilogx(w, np.angle(fr)*180/np.pi)
    plt.legend(('fem', 'fdt'))

plt.figure('Magnitude T/Pow x=-.01')
# mag, phase, omega = control.bode(control.frd(ft_pot(-.01,.0,0), w), w, Plot=False)
plt.semilogx(wf, Mf[0,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow x=+.01')
mag, phase, omega = control.bode(control.frd(ft_pot(.01,.0,0), w), w, Plot=False)
plt.semilogx(wf, Mf[1,:])
plt.semilogx(omega,20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow x=+.02')
mag, phase, omega = control.bode(control.frd(ft_pot(.02,.0,0), w), w, Plot=False)
plt.semilogx(wf, Mf[2,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow x=+.04')
mag, phase, omega = control.bode(control.frd(ft_pot(.04,.0,0), w), w, Plot=False)
plt.semilogx(wf, Mf[3,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow y=-.005')
mag, phase, omega = control.bode(control.frd(ft_pot(.0,-.005,0), w), w, Plot=False)
plt.semilogx(wf, Mf[4,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow y=.005')
mag, phase, omega = control.bode(control.frd(ft_pot(.0,.005,0), w), w, Plot=False)
plt.semilogx(wf, Mf[5,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow y=.01')
mag, phase, omega = control.bode(control.frd(ft_pot(.0,.01,0), w), w, Plot=False)
plt.semilogx(wf, Mf[6,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))
plt.figure('Magnitude T/Pow z=.005')
mag, phase, omega = control.bode(control.frd(ft_pot(.0,.0,0.005), w), w, Plot=False)
plt.semilogx(wf, Mf[7,:])
plt.semilogx(omega, 20*np.log10(mag))
plt.legend(('fem','fdt'))


# Phase
# idx = round(np.logspace(0,4,100));
# plt.semilogx(wf, af[0,:])
# af2 = af(:,idx);
# for i=1:4
#     for k=2:length(af2)
#         while af2(i,k)>(af2(i,k-1)+10)
#             af2(i,k)=af2(i,k)-360;
#         end
#     end
# end

X=np.load('output_data_vel.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi

Mf=None
af=None
U = np.fft.fft(X[:,2])
for i in range(3,X.shape[1]):
    H1 = U * np.fft.fft(X[:,i]) / U**2
    H1 = H1[:N]
    Mf = np.vstack((Mf, 20*np.log10(abs(H1)))) if Mf is not None else 20*np.log10(abs(H1))
    af = np.vstack((af, np.angle(H1)*180/np.pi)) if af is not None else np.angle(H1)*180/np.pi


plt.figure('Magnitude T/Vel')
for i in range(X.shape[1]-3):
    plt.semilogx(wf, Mf[i,:])
plt.legend(('x=-.01','x=.01','x=.02','x=.04','y=-.01','y=.005','y=.01','z=.005'))