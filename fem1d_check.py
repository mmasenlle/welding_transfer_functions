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
Tinf = 0
# piece dimensions
Lx=.4
# input source point
xi=.05
# output point from input source point
x1,x2=.01,.02
# inputs
q, v = 2e7, 0.002

# spatial discretization
nx=101
X=np.linspace(0,Lx,nx)

ops=np.load('output_points_1d_vel.npy')
x1,x2=ops[0],ops[1]

xo1,xo2=xi+x1,xi+x2

# assemble matrices
K=np.zeros((nx,nx))
C=np.zeros((nx,nx))
f=np.zeros(nx)
vu=np.array([[1],[1]])
NN2=np.array([[2,1],[1,2]])/6
lx=X[1]-X[0]
no1=None
no2=None

for i in range(nx-1):
    nds=np.array([i,i+1])
    AA=np.array([[1,X[i]],[1,X[i+1]]])
    abcd = np.linalg.solve(AA, np.eye(2))
    B = np.array([[abcd[1,0]], [abcd[1,1]]])
    Kek = (B @ B.T) * k * lx * A
    Kev = (vu @ B.T) * v * lx * A * cp * rho/2
    Keh = h * P * lx * NN2

    K[np.ix_(nds,nds)] += Kek + Kev + Keh
    C[np.ix_(nds,nds)] += np.eye(2) * cp * lx * A * rho/2
    f[nds] += h * P * Tinf * lx * np.ones(2) / 2

    if xi >= X[i] and xi < X[i+1]:
        f[i] += (abcd[0,0] + abcd[1,0]*xi) * q
        f[i+1] += (abcd[0, 1] + abcd[1, 1] * xi) * q
    if xo1 >= X[i] and xo1 < X[i+1]:
        no1 = np.array([[i,(X[i+1]-xo1)/lx],[i+1,(xo1-X[i])/lx]])
    if xo2 >= X[i] and xo2 < X[i+1]:
        no2 = np.array([[i,(X[i+1]-xo2)/lx],[i+1,(xo2-X[i])/lx]])

# force dirichlet conditions
bcwt = np.trace(K)/nx
K[0,:] = 0
K[:,0] = 0
K[0,0] = bcwt


# solve stationary system
T = np.linalg.solve(K, f)
if False:
    # plot temperatures on the source path
    plt.figure('Temperature profile')
    plt.plot(X,T)

if False:
    # assembled system from fenics
    A1=np.load('A1.npy')
    b1=np.load('b1.npy')
    u1=np.load('u1.npy')
    plt.plot(Lx-X,u1)
    # solved here
    TT = np.linalg.solve(A1, b1)
    plt.plot(Lx-X,TT)


################ ft
w = np.logspace(-4,1,200)
s = 1j*w
b=P*h*a/(k*A)
ft_pot = lambda x: (1/(rho*cp)*( np.exp((1/(2*a)*(x*v-np.sqrt(x**2*(v**2 + 4*a*b + 4*a*s))))))/np.sqrt(v**2 + 4*a*b + 4*a*s))
# control.bode_plot([control.frd(ft_pot(x1), w),control.frd(ft_pot(x2), w)], w, dB=True)


Asys=np.linalg.solve(C, -K)
Bsys=np.array([np.linalg.solve(C, f/q)]).T
Csys1=Csys2=np.zeros(nx)
Csys1[int(no1[0,0])],Csys1[int(no1[1,0])] = no1[0,1],no1[1,1]
Csys2[int(no2[0,0])],Csys2[int(no2[1,0])] = no2[0,1],no2[1,1]
Dsys=0
qsys1,qsys2 = control.ss(Asys,Bsys,Csys1,Dsys),control.ss(Asys,Bsys,Csys2,Dsys)
if False:
    plt.figure('Bode pot 1')
    control.bode_plot([qsys1,control.frd(ft_pot(x1),w)], w, dB=True)
    plt.legend(['qsys1','ft1'])
    plt.figure('Bode pot 2')
    control.bode_plot([qsys2,control.frd(ft_pot(x2),w)], w, dB=True)
    plt.legend(['qsys2','ft2'])


ft_vel = lambda x: (q/k*(np.exp((v*x)/(2*a))*(np.exp(-(x*np.sqrt(v**2 + 4*a*b))/(2*a))\
    - np.exp(-(x*np.sqrt(v**2 + 4*a*b + 4*a*s))/(2*a)) - (v*np.exp(-(x*np.sqrt(v**2 + 4*a*b))/(2*a)))\
    /np.sqrt(v**2 + 4*a*b) + (v*np.exp(-(x*np.sqrt(v**2 + 4*a*b + 4*a*s))/(2*a)))\
    /np.sqrt(v**2 + 4*a*b + 4*a*s)))/(2*s))
# control.bode_plot([control.frd(ft_vel(x1), w),control.frd(ft_vel(x2), w)], w, dB=True)

dT = np.concatenate((np.diff(T),T[-1:]))
Bvel=np.array([np.linalg.solve(C, dT)]).T*cp*rho/2
qvel1,qvel2 = control.ss(Asys,Bvel,Csys1,Dsys),control.ss(Asys,Bvel,Csys2,Dsys)
if False:
    plt.figure('Bode vel 1')
    control.bode_plot([qvel1,control.frd(ft_vel(x1),w)], w, dB=True)
    plt.legend(['qvel1','ft1'])
    plt.figure('Bode vel 2')
    control.bode_plot([qvel2,control.frd(ft_vel(x2),w)], w, dB=True)
    plt.legend(['qvel2','ft2'])

## data to bode

X=np.load('output_data_pow_1d.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,int(np.log10(N)),100).astype(int)
w = wf[idx]
s = 1j*w

ft_pot = lambda x: (1/(rho*cp)*( np.exp((1/(2*a)*(x*v-np.sqrt(x**2*(v**2 + 4*a*b + 4*a*s))))))/np.sqrt(v**2 + 4*a*b + 4*a*s))

U = np.fft.fft(X[:,2])
# for i in range(3,X.shape[1]):
#     H1 = U * np.fft.fft(X[:,i]) / U**2
#     plt.figure('Bode T/Pow ' + str(ops[i - 3]))
#     control.bode_plot((control.frd(H1[idx], w),control.frd(ft_pot(*ops[i - 3]), w)),w, dB=True)
#     plt.legend(('fem', 'fdt'))
H1 = U * np.fft.fft(X[:,3]) / U**2
H2 = U * np.fft.fft(X[:,4]) / U**2

plt.figure('Bode vel 1')
control.bode_plot([qsys1, control.frd(ft_pot(x1), w), control.frd(H1[idx], w)], w, dB=True)
plt.legend(['qsys1', 'ft1', 'fem1'])
plt.figure('Bode vel 2')
control.bode_plot([qsys2, control.frd(ft_pot(x2), w), control.frd(H2[idx], w)], w, dB=True)
plt.legend(['qsys2', 'ft2', 'fem2'])


# ops=np.load('output_points_1d_vel.npy')

X=np.load('output_data_vel_1d.npy')

# plt.plot(X[:1000,2])
dt = X[1,0]-X[0,0]
N = int(X.shape[0]/2)
freqs = np.arange(N)/(2*N*dt)
wf = freqs*2*np.pi
idx = np.logspace(0,int(np.log10(N)),100).astype(int)
w = wf[idx]
s = 1j*w

ft_vel = lambda x: (q/k*(np.exp((v*x)/(2*a))*(np.exp(-(x*np.sqrt(v**2 + 4*a*b))/(2*a))\
    - np.exp(-(x*np.sqrt(v**2 + 4*a*b + 4*a*s))/(2*a)) - (v*np.exp(-(x*np.sqrt(v**2 + 4*a*b))/(2*a)))\
    /np.sqrt(v**2 + 4*a*b) + (v*np.exp(-(x*np.sqrt(v**2 + 4*a*b + 4*a*s))/(2*a)))\
    /np.sqrt(v**2 + 4*a*b + 4*a*s)))/(2*s))

U = np.fft.fft(X[:,2])
# for i in range(3,X.shape[1]):
#     H1 = U * np.fft.fft(X[:,i]) / U**2
#     plt.figure('Bode T/Pow ' + str(ops[i - 3]))
#     control.bode_plot((control.frd(H1[idx], w),control.frd(ft_pot(*ops[i - 3]), w)),w, dB=True)
#     plt.legend(('fem', 'fdt'))
H1 = U * np.fft.fft(X[:,3]) / U**2
H2 = U * np.fft.fft(X[:,4]) / U**2

plt.figure('Bode vel 1')
control.bode_plot([qvel1, control.frd(ft_vel(x1), w), control.frd(H1[idx], w)], w, dB=True)
plt.legend(['qvel1', 'ft1', 'fem1'])
plt.figure('Bode vel 2')
control.bode_plot([qvel2, control.frd(ft_vel(x2), w), control.frd(H2[idx], w)], w, dB=True)
plt.legend(['qvel2', 'ft2', 'fem2'])