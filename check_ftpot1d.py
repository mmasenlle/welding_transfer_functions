import numpy as np
import matplotlib.pyplot as plt
import control,scipy

# material constants
k=24.0
rho=7925
cp=460
a=k/(rho*cp)
A=1
h=100000
P=2*np.pi*np.sqrt(A/np.pi)
b=P*h*a/(k*A)
Tinf = 0
# piece dimensions
Lx=.4
# input source point
xi=.05
# output point from input source point
# x1,x2=.01,.02
# inputs
q, v = 2e7, 0.002

# space state system from fem matrices (fenics1_test.py)
class Fem1d_fenics:
    def __init__(self):
        self.X = np.load('XX.npy')
        K = np.load('KK.npy')
        C = np.load('CC.npy')
        f = np.load('ff.npy')
        T0 = np.load('T0.npy')
        self.A = np.linalg.solve(C, -K)
        self.B = np.array([np.linalg.solve(C, f / q)]).T
        dT = np.concatenate((np.diff(T0), T0[-1:]))
        Bv = np.array([np.linalg.solve(C, dT)]).T * cp * rho / 2
    def get_ss(self, xo):
        C = np.zeros(self.X.size)
        C[np.where(self.X == xo)[0][0]] = 1
        return control.ss(self.A, self.B, C, 0)


fem1df = Fem1d_fenics()

w = np.logspace(-3, 3, 100)
s = 1j * w
ft_pot = lambda x: (1/(rho*cp)*( np.exp((1/(2*a)*(x*v-np.sqrt(x**2*(v**2 + 4*a*b + 4*a*s))))))/np.sqrt(v**2 + 4*a*b + 4*a*s))

for xo in (.015,.02,.03,.04):
    plt.figure('Bode pot xo=' + str(xo))
    control.bode_plot((control.frd(ft_pot(xo), w), fem1df.get_ss(xo + xi)), w, dB=True)
