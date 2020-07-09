# docker run -ti -v C:/src/sympy:/home/fenics/shared -w /home/fenics/shared fenics_ctrl
from dolfin import *
import numpy as np

# https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.html
# Create mesh and function space
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 32, 16, 4)
# mesh = Mesh('meshes/box2.xml')
# np.save('mesh_X',mesh.coordinates())
# np.save('mesh_cells',mesh.cells())


V = FunctionSpace(mesh, "CG", 2)
# V = FunctionSpace(mesh, "Q", 1)
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

k = Constant(24)
# k = Constant(30)
rho = 7925
Cp = 460

vel = Constant((0.005,0,0))

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

Pot=1500
delta = PointSource(V, Point(.02,.02,.0), Pot)

K, f = assemble_system(a1, L1, bc)

T0 = Function(V)
delta.apply(f)
solve(K, T0.vector(), f)

a = rho*Cp*u*v*dx
C = assemble(a)
# bc.apply(C)

# gT0 = project(T0.dx(0), V)

a = rho*Cp*inner(Constant((1,0,0)), nabla_grad(u))*v*dx
Kv = assemble(a)
# bc.apply(Kv)

print('About to solve AA')
AA3d = - np.linalg.solve(C.array(), K.array())
BB3d = np.array([np.linalg.solve(C.array(), f.get_local() / Pot)]).T
BB3dv = - np.array([np.linalg.solve(C.array(), Kv.array() @ T0.vector().get_local())]).T

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

print('About to store XX, AA, ...')
np.save('XX3d', dof_coordinates)
np.save('AA3d', AA3d)
np.save('BB3d', BB3d)
np.save('BB3dv', BB3dv)
np.save('T03d', np.array([T0.vector().get_local()]).T)

print('About to create T03d.pvd')
ofile = File("T03d.pvd")
ofile << T0
