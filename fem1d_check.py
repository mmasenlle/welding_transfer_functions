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

ops=np.load('output_points_1d_vel.npy')
x1,x2=ops[0],ops[1]

xo1,xo2=xi+x1,xi+x2

class Fem1d:
    def __init__(self, nx = 81):
        # spatial discretization
        X = self.X = np.linspace(0, Lx, nx)
        # assemble matrices
        self.K=np.zeros((nx,nx))
        self.C=np.zeros((nx,nx))
        self.f=np.zeros(nx)
        vu=np.array([[1],[1]])
        NN2=np.array([[2,1],[1,2]])/6
        lx=X[1]-X[0]
        self.no1=None
        self.no2=None

        for i in range(nx-1):
            nds=np.array([i,i+1])
            AA=np.array([[1,X[i]],[1,X[i+1]]])
            abcd = np.linalg.solve(AA, np.eye(2))
            B = np.array([[abcd[1,0]], [abcd[1,1]]])
            Kek = (B @ B.T) * k * lx * A
            Kev = (vu @ B.T) * v * lx * A * cp * rho/2
            Keh = h * P * lx * NN2

            self.K[np.ix_(nds,nds)] += Kek + Kev + Keh
            self.C[np.ix_(nds,nds)] += NN2 * cp * lx * A * rho
            self.f[nds] += h * P * Tinf * lx * np.ones(2) / 2

            if xi >= X[i] and xi < X[i+1]:
                self.f[i] += (abcd[0,0] + abcd[1,0]*xi) * q
                self.f[i+1] += (abcd[0, 1] + abcd[1, 1] * xi) * q
            if xo1 >= X[i] and xo1 < X[i+1]:
                self.no1 = np.array([[i,(X[i+1]-xo1)/lx],[i+1,(xo1-X[i])/lx]])
            if xo2 >= X[i] and xo2 < X[i+1]:
                self.no2 = np.array([[i,(X[i+1]-xo2)/lx],[i+1,(xo2-X[i])/lx]])

        # force dirichlet conditions
        bcwt = np.trace(self.K)/nx
        self.K[0,:] = 0
        self.K[:,0] = 0
        self.K[0,0] = 1

        # solve stationary system
        self.T0 = np.linalg.solve(self.K, self.f)

    def plot(self):
        plt.figure('Fem1d T0')
        plt.plot(self.X, self.T0)

    def init_ss(self):
        A = np.linalg.solve(self.C, -self.K)
        B = np.array([np.linalg.solve(self.C, self.f / q)]).T
        C1 = np.zeros(self.X.size)
        C2 = np.zeros(self.X.size)
        no1,no2 = self.no1,self.no2
        C1[int(no1[0, 0])], C1[int(no1[1, 0])] = no1[0, 1], no1[1, 1]
        C2[int(no2[0, 0])], C2[int(no2[1, 0])] = no2[0, 1], no2[1, 1]
        D = 0
        self.ss_pot1, self.ss_pot2 = control.ss(A, B, C1, D), control.ss(A, B, C2, D)
        dT = np.concatenate((np.diff(self.T0), self.T0[-1:]))
        Bvel = np.array([np.linalg.solve(self.C, dT)]).T * cp * rho
        self.ss_vel1, self.ss_vel2 = control.ss(A, Bvel, C1, D), control.ss(A, Bvel, C2, D)


# space state system from fem matrices (fenics1_test.py)
class Fem1d_fenics:
    def __init__(self):
        self.X = np.load('XX.npy')
        self.K = np.load('KK.npy')
        self.C = np.load('CC.npy')
        self.f = np.load('ff.npy')
        self.T0 = np.load('T0.npy')
    def init_ss(self):
        A = np.linalg.solve(self.C, -self.K)
        B = np.array([np.linalg.solve(self.C, self.f / q)]).T
        C1 = np.zeros(self.X.size)
        C2 = np.zeros(self.X.size)
        C1[np.where(self.X == xo1)[0][0]] = 1
        C2[np.where(self.X == xo2)[0][0]] = 1
        D = 0
        self.ss_pot1, self.ss_pot2 = control.ss(A, B, C1, D), control.ss(A, B, C2, D)
        dT = np.concatenate((np.diff(self.T0), self.T0[-1:]))
        Bvel = np.array([np.linalg.solve(self.C, dT)]).T * cp * rho
        self.ss_vel1, self.ss_vel2 = control.ss(A, Bvel, C1, D), control.ss(A, Bvel, C2, D)


fem1dp = Fem1d()
fem1df = Fem1d_fenics()

if False:
    plt.figure('Fem1d T0')
    plt.plot(fem1dp.X, fem1dp.T0)
    plt.plot(fem1df.X, fem1df.T0)
    plt.legend(('fem1dp','fem1df'))


## Frequency response from simulation data
class FrdSimulPow:
    def __init__(self):
        X=np.load('output_data_pow_1d.npy')

        # plt.plot(X[:1000,2])
        dt = X[1,0]-X[0,0]
        N = int(X.shape[0]/2)
        freqs = np.arange(N)/(2*N*dt)
        wf = freqs*2*np.pi
        idx = np.logspace(0,int(np.log10(N))-1,100).astype(int)
        self.w = wf[idx]

        U = np.fft.fft(X[:,1])
        H1 = U * np.fft.fft(X[:,3]) / U**2
        H2 = U * np.fft.fft(X[:,4]) / U**2
        self.H1 = H1[idx]
        self.H2 = H2[idx]

class FrdSimulVel:
    def __init__(self):
        X = np.load('output_data_vel_1d.npy')
        dt = X[1, 0] - X[0, 0]
        N = int(X.shape[0] / 2)
        freqs = np.arange(N) / (2 * N * dt)
        wf = freqs * 2 * np.pi
        idx = np.logspace(0, int(np.log10(N)) - 1, 100).astype(int)
        self.w = wf[idx]

        U = np.fft.fft(X[:, 2])
        H1 = U * np.fft.fft(X[:, 3]) / U ** 2
        H2 = U * np.fft.fft(X[:, 4]) / U ** 2
        self.H1 = H1[idx]
        self.H2 = H2[idx]

fem1dp.init_ss()
fem1df.init_ss()

frdPow = FrdSimulPow()
w = frdPow.w
s = 1j * w
ft_pot = lambda x: (1/(rho*cp)*( np.exp((1/(2*a)*(x*v-np.sqrt(x**2*(v**2 + 4*a*b + 4*a*s))))))/np.sqrt(v**2 + 4*a*b + 4*a*s))
plt.figure('Bode pot 1')
control.bode_plot([control.frd(frdPow.H1, w), control.frd(ft_pot(x1), w), fem1dp.ss_pot1, fem1df.ss_pot1], w, dB=True)
plt.legend(['fem', 'ft', 'ssp', 'ssf'])
plt.figure('Bode pot 2')
control.bode_plot([control.frd(frdPow.H2, w), control.frd(ft_pot(x2), w), fem1dp.ss_pot2, fem1df.ss_pot2], w, dB=True)
plt.legend(['fem', 'ft', 'ssp', 'ssf'])


frdVel = FrdSimulVel()
w = frdVel.w
s = 1j * w
ft_vel = lambda x: (q / k * (np.exp((v * x) / (2 * a)) * (np.exp(-(x * np.sqrt(v ** 2 + 4 * a * b)) / (2 * a)) - np.exp(
                    -(x * np.sqrt(v ** 2 + 4 * a * b + 4 * a * s)) / (2 * a)) - (v * np.exp(
                    -(x * np.sqrt(v ** 2 + 4 * a * b)) / (2 * a))) / np.sqrt(v ** 2 + 4 * a * b) + (v * np.exp(
                    -(x * np.sqrt(v ** 2 + 4 * a * b + 4 * a * s)) / (2 * a)))/ np.sqrt(v ** 2 + 4 * a * b + 4 * a * s))) / (2 * s))
plt.figure('Bode vel 1')
control.bode_plot([control.frd(frdVel.H1, w), control.frd(ft_vel(x1), w), fem1dp.ss_vel1, fem1df.ss_vel1], w, dB=True)
plt.legend(['fem', 'ft', 'ssp', 'ssf'])
plt.figure('Bode vel 2')
control.bode_plot([control.frd(frdVel.H2, w), control.frd(ft_vel(x2), w), fem1dp.ss_vel2, fem1df.ss_vel2], w, dB=True)
plt.legend(['fem', 'ft', 'ssp', 'ssf'])

