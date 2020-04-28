#  docker run -ti -v C:/src/sympy:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current

from dolfin import *
import numpy as np

# Create mesh and function space
mesh = IntervalMesh(80,0,0.4)
V = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, DirichletBoundary())


Pot=2e7
k = Constant(24)
rho = Constant(7925)
#Pot=2
#k = Constant(0.000024)
#rho = Constant(0.007925)
Cp = Constant(460)

h = Constant(100000)
#h = Constant(0.001)
Tinf = Constant(0)
P = Constant(3.5449)

v_nom = 0.002
vel = Constant(v_nom)

ipx = .05

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*vel*u.dx(0)*v*dx + P*h*inner(u,v)*dx
L1 = Constant(0)*v*dx + P*h*Tinf*v*dx

A1, b1 = assemble_system(a1, L1, bc)
delta = PointSource(V, Point(.05,), Pot)
delta.apply(b1)

u1 = Function(V)

# estacionario
solve(A1, u1.vector(), b1)
print ("max: ", max(u1.vector()),"(K), T(x=0.2): ", u1.vector()[50]) # 137.15

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

# check steady state temperature field
# np.save('X_2d', dof_coordinates)
# np.save('T0_2d', u1.vector().get_local())
# exit(0)
output_points0 = (.015,.02,.03,.04)
np.save('output_points_1d_pow', output_points0)
output_points = np.array(output_points0) + np.array(ipx)


opidx=[]
for op in output_points:
    opidx.append(np.argwhere(np.all((dof_coordinates-np.array(op))==0, axis=1))[0][0])
    print("T(", dof_coordinates[opidx[-1]], ",t=0):", u1.vector()[opidx[-1]], "K")
print("Tmax:", np.max(u1.vector().get_local()), "K")

# transitorio
u_prev = u1 #interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

b = None
uf = Function(V)

t = 0.0
# print ("fenics_bode_inc=[", file=sys.stderr)
# print (repr(t), repr(Pot), repr(v_nom), repr(u1.vector()[idxT1]), repr(u1.vector()[idxT2]), repr(u1.vector()[idxT3]), repr(u1.vector()[idxT4]), file=sys.stderr)
output_data = np.array([t,Pot,v_nom]+[u1.vector()[i] for i in opidx])

N = 1024*256
dt = .0002
T = N*dt
freq = .0
df = .0002
for ii in range(N-1):
    t += dt
    #vt = v_nom +(v_nom*.1*np.sin(2*pi*freq*t))
    vel = Constant(v_nom)
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*vel*u.dx(0)*v*dx + dt*P*h*u*v*dx
    L = rho*Cp*u_prev*v*dx + dt*P*h*Tinf*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    bc.apply(A, b)
    Pt = Pot +(Pot * np.sin(2*pi*freq*t))
    delta = PointSource(V, Point(.05,), dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    # print (repr(t), repr(Pt), repr(Pot), repr(uf.vector()[idxT1]), repr(uf.vector()[idxT2]), repr(uf.vector()[idxT3]), repr(uf.vector()[idxT4]), file=sys.stderr)
    output_data = np.vstack((output_data,[t, Pt, v_nom]+[uf.vector()[i] for i in opidx]))
    u_prev.assign(uf)
    freq += df
    if ii % 77 == 0:
        print("t=", round(t,3), 'of', T, 'seconds')

# print ("];", file=sys.stderr)
np.save('output_data_pow_1d', output_data)
for opi in opidx:
    print("T(", dof_coordinates[opi], ",t=",T,"):", uf.vector()[opi], "K")
print("Tmax:", np.max(uf.vector().get_local()), "K")