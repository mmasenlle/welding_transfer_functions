from dolfin import *

# Create mesh and function space
mesh = IntervalMesh(100,0,0.4)
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

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*vel*u.dx(0)*v*dx + P*h*inner(u,v)*dx
L1 = Constant(0)*v*dx + P*h*Tinf*v*dx

A1, b1 = assemble_system(a1, L1, bc)
delta = PointSource(V, Point(0.05,), Pot)
delta.apply(b1)

u1 = Function(V)
# estacionario
solve(A1, u1.vector(), b1)

import numpy as np
np.save('A1', A1.array())
np.save('b1', b1.get_local())
np.save('u1', u1.vector().get_local())