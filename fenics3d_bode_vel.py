#  docker run -ti -v C:/src/sympy:/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current
import sys
from dolfin import *
import numpy as np
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

output_points_vel = ((-.01,0,0),(-.005,0,0),(.005,0,0),(.01,0,0),(.015,0,0),(.02,0,0),(.03,0,0),(.04,0,0),
                 (0,-.003,0),(0,.003,0),(0,-.005,0),(0,.005,0),(0,-.007,0),(0,.007,0))
                 # (.01,.005,0),(-.01,.005,0),(.01,.01,0),(.01,0,.005))
# output_points_vel = ((.01,0,0),)
output_points_vel = ()
save_data = False
if MPI.COMM_WORLD.Get_rank() == 0 and save_data:
    np.save('output_points_vel_all', output_points_vel)
# output_points = np.array(output_points_vel) + np.array((ipx,ipy,0)) if len(output_points_vel) > 0 else []

bbt = mesh.bounding_box_tree()
o_points=[]
for o in output_points_vel:
    op = np.array(o) + np.array((ipx,ipy,0))
    if bbt.compute_first_entity_collision(Point(*op)) < mesh.num_cells():
        print(MPI.COMM_WORLD.Get_rank(), "T(", *op, ",t=0):", u1(*op), "K")
        o_points.append(op)
if save_data: np.save('output_points_vel_' + str(MPI.COMM_WORLD.Get_rank()), o_points)
print('MPI Rank:', MPI.COMM_WORLD.Get_rank(), '/', MPI.COMM_WORLD.Get_size(), 'points:', len(o_points))

# transitorio
u_prev = u1 #interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

b = None
uf = Function(V)

t = 0.0
# print ("fenics_bode_inc=[", file=sys.stderr)
# print (repr(t), repr(Pot), repr(v_nom), repr(u1.vector()[idxT1]), repr(u1.vector()[idxT2]), repr(u1.vector()[idxT3]), repr(u1.vector()[idxT4]), file=sys.stderr)
# output_data = np.array([t,Pot,v_nom]+[u1.vector()[i] for i in opidx])
output_data = np.array([t,Pot,v_nom]+[u1(*op) for op in o_points])

N = 100 #1024*512
dt = .001
T = N*dt
freq = .0
df = .0005
for ii in range(N-1):
    t += dt
    vt = v_nom +(v_nom*.1*np.sin(2*pi*freq*t))
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*inner(Constant((vt,0,0)), nabla_grad(u))*v*dx
    L = rho*Cp*u_prev*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    bc.apply(A, b)
    Pt = Pot#+(Pot * np.sin(2*pi*freq*t))
    delta = PointSource(V, inp, dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    # print (repr(t), repr(Pt), repr(Pot), repr(uf.vector()[idxT1]), repr(uf.vector()[idxT2]), repr(uf.vector()[idxT3]), repr(uf.vector()[idxT4]), file=sys.stderr)
    # output_data = np.vstack((output_data,[t, Pt, vt]+[uf.vector()[i] for i in opidx]))
    output_data = np.vstack((output_data, [t, Pt, vt] + [uf(*op) for op in o_points]))
    u_prev.assign(uf)
    freq += df
    if ii % 77 == 0:
        print("t=", round(t,3), 'of', T, 'seconds')

# print ("];", file=sys.stderr)
if save_data: np.save('output_data_vel_' + str(MPI.COMM_WORLD.Get_rank()), output_data)
for op in o_points:
    print("T(", *op, ",t=",T,"):", uf(*op), "K")
print("Tmax:", np.max(uf.vector().get_local()), "K")
