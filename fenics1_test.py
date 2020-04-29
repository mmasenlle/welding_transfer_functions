from dolfin import *

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

Cp = Constant(460)

h = Constant(100000)

Tinf = Constant(0)
P = Constant(3.5449)

v_nom = 0.002
vel = Constant(v_nom)

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*vel*u.dx(0)*v*dx + P*h*inner(u,v)*dx
L1 = Constant(0)*v*dx + P*h*Tinf*v*dx

A1, b1 = assemble_system(a1, L1, bc)
delta = PointSource(V, Point(0.05,), Pot)
delta.apply(b1)

u1 = Function(V)

solve(A1, u1.vector(), b1)

a = rho*Cp*u*v*dx
A = assemble(a)

import numpy as np
np.save('KK', A1.array())
np.save('ff', b1.get_local())
np.save('T0', u1.vector().get_local())
np.save('CC', A.array())
np.save('XX', V.tabulate_dof_coordinates())
