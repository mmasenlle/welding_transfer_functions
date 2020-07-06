from dolfin import *
import numpy as np

# Create mesh and function space
mesh = RectangleMesh(Point(0.0, 0.0), Point(.16, .04), 64, 16)
# np.save('mesh_X',mesh.coordinates())
# np.save('mesh_cells',mesh.cells())


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

Lz=.01
Pot=1000/Lz
ipx=.02
ipy=.02
inp = Point(ipx,ipy)

k = Constant(24)
rho = Constant(7925)
#Pot=2
#k = Constant(0.000024)
#rho = Constant(0.007925)
Cp = Constant(460)
h = Constant(1000)

v_nom = 0.002
vel = Constant((v_nom,0))

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx + 2*h*inner(u,v)*dx
L1 = Constant(0)*v*dx

delta = PointSource(V, inp, Pot)

K, f = assemble_system(a1, L1, bc)

T0 = Function(V)
delta.apply(f)
solve(K, T0.vector(), f)

a = rho*Cp*u*v*dx
C = assemble(a)
# bc.apply(C)

# gT0 = project(T0.dx(0), V)

a = rho*Cp*inner(Constant((1,0)), nabla_grad(u))*v*dx
Kv = assemble(a)
# bc.apply(Kv)

AA2d = - np.linalg.solve(C.array(), K.array())
BB2d = np.array([np.linalg.solve(C.array(), f.get_local() / Pot)]).T
BB2dv = - np.array([np.linalg.solve(C.array(), Kv.array() @ T0.vector().get_local())]).T

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

np.save('XX2d', dof_coordinates)
np.save('AA2d', AA2d)
np.save('BB2d', BB2d)
np.save('BB2dv', BB2dv)