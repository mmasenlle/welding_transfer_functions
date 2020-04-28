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
b=0

ops=np.load('output_points_2d_vel.npy')

X=np.load('output_data_vel_2d.npy')

# plt.plot(X[:10000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,int(np.log10(N)),100).astype(int)
# idx = np.logspace(0,4,100).astype(int)
idx[0] = 0
w = wf[idx]
s = 1j*w
s[0] = 1j*1e-10
K=k
q=100000

R = lambda x,y: np.sqrt(x**2+y**2)
D = lambda s: np.sqrt(4*a*s + 4*a*b + v**2)

ft_vel1 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a)))
ft_vel2 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
    x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/(4*a) - x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/(4*a)/(4*a/(4*a*b + v**2)*s + 1))
# ft_vel3 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
#                                 + x*D(0)*D(0)*kv(0,R(x,y)*D(0)/(2*a))/(4*a) - x*D(0)**4*kv(0,R(x,y)*D(0)/(2*a))/(4*a)/D(s)**2)
ft_vel3 = lambda x,y: ft_vel1(x,y) + ft_vel2(x,y)
# x=-0.01 a1 = 2*a
# x=-0.01 y 0.01 a1 = 2.35*a
# x=0.02 a1 = 5.35*a
# x=0.04 a1 = 11.7*a
a1 = 3.037916666666667*a
# a1 = lambda x: np.abs(x)*235*a
a1 = lambda x: x*(4*a*b + v**2)/(4*a*v)*4*a
ft_vel4 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
    x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/a1(x) - x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/a1(x)/(a1(x)/(4*a*b + v**2)*s + 1))
ft_vel5 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y)

ft_vel1_ = lambda x,y: q/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a)))



aa = 1
bb = 1
# aa/bb
ft_vel6 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a)) + kv(0,R(x,y)*D(0)/(2*a))*bb/(aa) - kv(0,R(x,y)*D(0)/(2*a))*bb**2/aa/(aa*x*s/v + 1*bb))
ft_vel61 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a)))
ft_vel62 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(kv(0,R(x,y)*D(0)/(2*a))*bb/(aa) - kv(0,R(x,y)*D(0)/(2*a))*bb**2/aa/(aa*x*s/v + 1*bb))

aa = .8
ft_vel7 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a)) + kv(0,R(x,y)*D(0)/(2*a))/aa - kv(0,R(x,y)*D(0)/(2*a))/aa/np.sqrt(aa*2*x*s/v + 1))

aaa=.7
ft_vel8 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
                                + x*D(0)**3*kv(0,R(x,y)*D(0)/(2*a))/np.sqrt(4*a*b + v**2)/(2*a*aaa) - x*D(0)**3*kv(0,R(x,y)*D(0)/(2*a))/np.sqrt(aaa*4*a*s + 4*a*b + v**2)/(2*a*aaa))


aa = 1
bb = aa/1
# aa/bb
ft_vel9 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a)) + kv(0,R(x,y)*D(0)/(2*a))*bb/2/(aa) - kv(0,R(x,y)*D(0)/(2*a))*bb**2/2/aa/(aa*x*s/v + 1*bb)**2)

ft_vel10 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a)) + x*kv(0,R(x,y)*D(0)/(2*a))*s/(x*s/v + 1))


aa = 1
ft_vel11 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a)) + x*kv(0,R(x,y)*D(0)/(2*a))*aa - x*kv(0,R(x,y)*D(0)/(2*a))*aa*(s + 1))

# the best so far
# ft_vel12 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
#         v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
#         + x/R(x,y)**2*2*a*kv(0,R(x,y)*D(s)/(2*a)) - x/R(x,y)**2*2*a*kv(0,R(x,y)*D(0)/(2*a))\
#         - x*D(0)/R(x,y)*kv(1,R(x,y)*D(s)/(2*a)) + x*D(0)/R(x,y)*kv(1,R(x,y)*D(0)/(2*a)))
ft_vel12 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        (v + x*2*a/R(x,y)**2)*(kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a)))\
        - x*D(0)/R(x,y)*(kv(1,R(x,y)*D(s)/(2*a)) - kv(1,R(x,y)*D(0)/(2*a))))
a1=lambda x,y: (v + x*2*a/R(x,y)**2)
a2=lambda x,y: x*D(0)/R(x,y)

aa=1.0
bb=.54/aa
ft_vel12a = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        (v + aa*bb*x*2*a/R(x,y)**2)*(kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a))))
ft_vel12b = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        - aa*x*D(0)/R(x,y)*(kv(1,R(x,y)*D(s)/(2*a)) - kv(1,R(x,y)*D(0)/(2*a))))
ft_vel12ab = lambda x,y: -np.sign(x)*(ft_vel12a(x,y)+ft_vel12b(x,y))

ft_vel13 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
        + x*D(0)**2*kv(0,R(x,y)*D(s)/(2*a))/(2*a) - x*D(0)**2*kv(0,R(x,y)*D(0)/(2*a))/(2*a)\
        - x*D(0)**3*kv(0,R(x,y)*D(s)/(2*a))/D(s)/(2*a) + x*D(0)**3*kv(0,R(x,y)*D(0)/(2*a))/D(0)/(2*a))

ft_vel14 = lambda x,y: q*v*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(
        kv(0,R(x,y)*D(s)/(2*a)) - kv(0,R(x,y)*D(0)/(2*a))/(x*s/v + 1))

# aa1=np.array((2,2.1,2.2,2.3,2.4,2.5,2.6))*a
# ft_vel_a1 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y,aa1[1])
# ft_vel_a2 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y,aa1[2])
# ft_vel_a3 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y,aa1[3])
# ft_vel_a4 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y,aa1[4])
# ft_vel_a5 = lambda x,y: ft_vel1(x,y) + ft_vel4(x,y,aa1[5])
# pole = (4*a*b + v**2)/a1
# pole_log = np.log10(pole)
# print('pole', pole, pole_log)
# ft_vel6 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
#                                 + x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/(4*a) - x*(4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/(4*a)/(4*a/(4*a*b + v**2)*s + 1))
# ft_vel7 = lambda x,y: q*np.exp(v*x/(2*a))/(4*np.pi*K*a*s)*(v*kv(0,R(x,y)*D(s)/(2*a)) - v*kv(0,R(x,y)*D(0)/(2*a))\
#                                 + (4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a)) - (4*a*b + v**2)*kv(0,R(x,y)*D(0)/(2*a))/(x**2/(4*a*b + v**2)*s + 1))


ft_vel0 = lambda x,y: q*(-v*R(x,y)*kv(1, R(x,y)*D(0)/(2*a)) + x*D(0)*kv(0, R(x,y)*D(0)/(2*a)))*np.exp(v*x/(2*a))/(4*np.pi*K*a*D(0))

if False:
    H1 = U * np.fft.fft(X[:, 4]) / U ** 2
    control.bode_plot((control.frd(H1[idx], w),control.frd(ft_vel(*ops[1]), w)),w, dB=True)
    ft_vel1 = lambda x, y: q * np.exp(v * x / (2 * a)) / (4 * np.pi * K * a * s) * (
                v * kv(0, R(x, y) * D(s) / (2 * a)) - v * kv(0, R(x, y) * D(0) / (2 * a))\
                + x * kv(0, R(x, y) * D(0) / (2 * a))/(R(x, y) * D(s)))
    control.bode_plot(control.frd(ft_vel1(.01,0), w), w, dB=True)

U = np.fft.fft(X[:,2])
for i in range(3,X.shape[1]):
    H1 = U * np.fft.fft(X[:,i]) / U**2
    plt.figure('Bode T/Pow ' + str(ops[i - 3]))
    control.bode_plot((control.frd(H1[idx], w), control.frd(ft_vel12(*ops[i - 3]), w)),w, dB=True)
                       # control.frd(ft_vel_a1(*ops[i - 3]), w),control.frd(ft_vel_a2(*ops[i - 3]), w),
                       # control.frd(ft_vel3(*ops[i - 3]), w), control.frd(ft_vel6(*ops[i - 3]), w),
                       # control.frd(ft_vel13(*ops[i - 3]), w),
                       # control.frd(ft_vel12(*ops[i - 3]), w), control.frd(ft_vel12ab(*ops[i - 3]), w)),w, dB=True)
    plt.legend(('fem', 'fdt12', 'fdt12ab'))
    plt.axvline(v/ops[i-3][0], color='yellow', linestyle='dashed')
    print(ft_vel0(*ops[i - 3]), np.abs(ft_vel1_(*ops[i - 3])[1]), np.abs(H1[1]))

if False:
    x1=0.01
    y1=.08
    control.bode_plot((control.frd(ft_vel61(x1, y1), w),control.frd(ft_vel62(x1, y1), w),control.frd(ft_vel61(x1, y1)+ft_vel62(x1, y1), w)), w, dB=True)
    plt.legend(('1', '2', '1+2'))
