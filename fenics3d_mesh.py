from dolfin import *
import numpy as np

# Create mesh and function space
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 32, 16, 4)
np.save('mesh_X',mesh.coordinates())
np.save('mesh_cells',mesh.cells())


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

k = Constant(24)
rho = Constant(7925)
Cp = Constant(460)

vel = Constant((0.002,0,0))

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

Pot=100
delta = PointSource(V, Point(.02,.02,.0), Pot)

A1, b1 = assemble_system(a1, L1, bc)

u1 = Function(V)
delta.apply(b1)
solve(A1, u1.vector(), b1)

np.save('A13d',A1.array())
np.save('b13d',b1.get_local())
np.save('u13d',u1.vector().get_local())

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))
np.save('dof_coordinates', dof_coordinates)