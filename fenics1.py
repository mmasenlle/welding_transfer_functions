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
print ("max: ", max(u1.vector()),"(K), T(x=0.2): ", u1.vector()[50]) # 137.15

# Control
idxT1 = 100 - 15 # from Matlab
refT1 = u1.vector()[idxT1]
#idxT2 = 1000 - 375 # from Matlab
idxT2 = 100 - 17 # from Matlab
refT2 = u1.vector()[idxT2]
idxT3 = 100 - 20 # from Matlab
refT3 = u1.vector()[idxT3]
idxT4 = 100 - 22 # from Matlab
refT4 = u1.vector()[idxT4]

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))
dof_x = dof_coordinates[:, 0]
print ("T(x=",dof_x[idxT1],",t=0): ", refT1, " K")
print ("T(x=",dof_x[idxT2],",t=0): ", refT2, " K")
print ("T(x=",dof_x[idxT3],",t=0): ", refT3, " K")
print ("T(x=",dof_x[idxT4],",t=0): ", refT4, " K")