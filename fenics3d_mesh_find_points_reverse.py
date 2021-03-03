# docker run -ti -v C:/src/sympy:/home/fenics/shared fenics1
from dolfin import *
import numpy as np

mesh = Mesh('meshes/box_reverse.xml')
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

# Thermal analysis of the arc welding process: Part I. General solutions, Table 1
k = Constant(60)
rho = 7870
Cp = 2719

vel = Constant((0.005,0,0))

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

Pot=7500
delta = PointSource(V, Point(.06,.02,.0), Pot)

K, f = assemble_system(a1, L1, bc)

T0 = Function(V)
delta.apply(f)
solve(K, T0.vector(), f)

points=[]
isoterms = (1200,) #1000,800,600)
p1=0.06
p2=0.02
for isot in isoterms:
    while T0(Point(p1, 0.02, 0)) > isot:
        # p1 += .000001
        p1 -= .000001
    points.append(p1)
    print('T(%f,0.02,0): %f %f' % (p1, T0(Point(p1, 0.02, 0)), 0.06-p1))
    while T0(Point(0.06, p2, 0)) > isot:
        p2 += .000001
    points.append(p2)
    print('T(0.06,%f,0): %f %f' % (p2, T0(Point(0.06, p2, 0)), p2-0.02))
    print(str((p2-0.02)/(0.06-p1)))

# print('points:',points)