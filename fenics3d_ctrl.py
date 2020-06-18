#  docker run -ti -v C:/src/sympy:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current
import sys
from dolfin import *
import numpy as np
import control
import filter1
from mpi4py import MPI

# Create mesh and function space
# mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 320, 160, 40)
mesh = Mesh('meshes/box.xml')
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

Pot=1500
ipx=.02
ipy=.02
inp = Point(ipx,ipy,.0)

k = Constant(24)
rho = Constant(7925)
#Pot=2
#k = Constant(0.000024)
#rho = Constant(0.007925)
Cp = Constant(460)

v_nom = 0.005
vel = Constant((v_nom,0,0))

# a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*vel*u.dx(0)*v*dx + P*h*inner(u,v)*dx
# L1 = Constant(0)*v*dx + P*h*Tinf*v*dx
a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

A1, b1 = assemble_system(a1, L1, bc)
delta = PointSource(V, inp, Pot)
delta.apply(b1)

u1 = Function(V)

# estacionario
solve(A1, u1.vector(), b1)
# print ("max: ", max(u1.vector()),"(K), T(x=0.2): ", u1.vector()[500]) # 137.15

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

T1ref=1000
T2ref=400
p1=Point(0.03,0.02,0)
p2=Point(0.02,0.025,0)

# transitorio
u_prev = u1 #interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

b = None
uf = Function(V)



t = 0.0
dt = .001
T = 50

# Controllers
ctrl1_d = control.sample_system(control.tf(1, (1, 0)), dt)
ctrl1 = filter1.Filter1((0,ctrl1_d.num[0][0][0]),ctrl1_d.den[0][0])
ctrl2_d = control.sample_system(control.tf(-0.0002, (1, 0)), dt)
ctrl2 = filter1.Filter1((0,ctrl2_d.num[0][0][0]),ctrl2_d.den[0][0])

output_data = np.array([t,Pot,v_nom,u1(p1),u1(p2)])

tn=1
T1,T2=u_prev(p1),u_prev(p2)
stable=0
while t <= T:
    t += dt
    Pt = Pot + ctrl1.step(T1ref - T1)
    vt = v_nom + ctrl2.step(T2ref - T2)
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*inner(Constant((vt,0,0)), nabla_grad(u))*v*dx
    L = rho*Cp*u_prev*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    bc.apply(A, b)
    delta = PointSource(V, inp, dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    u_prev.assign(uf)
    T1, T2 = u_prev(p1), u_prev(p2)
    output_data = np.vstack((output_data, [t, Pt, vt, uf(p1), uf(p2)]))
    if T1**2+T2**2 < .01:
        stable += 1
        if stable > 100:
            print("Output stable t=", round(t, 3), 'of', T, 'seconds')
            break
    else:
        stable = 0
    if t >= tn:
        print("t=", round(t,3), 'of', T, 'seconds')
        tn += 1

# print ("];", file=sys.stderr)
# np.save('output_data_vel_' + str(MPI.COMM_WORLD.Get_rank()), output_data)
# for op in o_points:
#     print("T(", *op, ",t=",T,"):", uf(*op), "K")
np.save('ctrl_data', output_data)
print("Tmax:", np.max(uf.vector().get_local()), "K")
