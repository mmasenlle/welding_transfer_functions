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

XX = np.load('XX.npy')
KK = np.load('KK.npy')
CC = np.load('CC.npy')
ff = np.load('ff.npy')
T0 = np.load('T0.npy')

AA = np.linalg.solve(CC, -KK)
BB = np.array([np.linalg.solve(CC, ff / q)]).T

nx = XX.size
X = np.linspace(0, Lx, nx)
# assemble matrices
K=np.zeros((nx,nx))
C=np.zeros((nx,nx))
f=np.zeros(nx)
vu=np.array([[1],[1]])
NN2=np.array([[2,1],[1,2]])/6
lx=X[1]-X[0]

for i in range(nx-1):
    nds=np.array([i,i+1])
    AA=np.array([[1,X[i]],[1,X[i+1]]])
    abcd = np.linalg.solve(AA, np.eye(2))
    B = np.array([[abcd[1,0]], [abcd[1,1]]])
    Kek = (B @ B.T) * k * lx * A
    Kev = (vu @ B.T) * v * lx * A * cp * rho/2
    Keh = h * P * lx * NN2

    K[np.ix_(nds,nds)] += Kek + Kev + Keh
    C[np.ix_(nds,nds)] += NN2 * cp * lx * A * rho
    f[nds] += h * P * Tinf * lx * np.ones(2) / 2

    if xi >= X[i] and xi < X[i+1]:
        f[i] += (abcd[0,0] + abcd[1,0]*xi) * q
        f[i+1] += (abcd[0, 1] + abcd[1, 1] * xi) * q

# force dirichlet conditions
bcwt = np.trace(K)/nx
K[0,:] = 0
K[:,0] = 0
K[0,0] = 1

KK-np.flip(K)


