# docker run -ti -v C:/src/sympy:/home/fenics/shared -w /home/fenics/shared fenics1
from dolfin import *
import numpy as np

# Create mesh and function space
# mesh = RectangleMesh(Point(0.0, 0.0), Point(.16, .04), 64, 16)
# mesh = RectangleMesh(Point(0.0, 0.0), Point(.16, .04), 256, 64)
mesh = Mesh('meshes/box2d_2.xml')
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
Pot=2000/Lz
ipx=.02
ipy=.02
inp = Point(ipx,ipy)

k = Constant(24)
rho = Constant(7925)
#Pot=2
#k = Constant(0.000024)
#rho = Constant(0.007925)
Cp = Constant(460)
h = Constant(0) #Constant(1000)

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

# np.save('XX2d', dof_coordinates)
# np.save('AA2d', AA2d)
# np.save('BB2d', BB2d)
# np.save('BB2dv', BB2dv)


def get_C(xo):
    C = np.zeros(dof_coordinates.shape[0])
    X_xo = dof_coordinates - np.array(xo)
    idx = np.argwhere(np.all((X_xo) == 0, axis=1))
    if idx.size:
        C[idx[0][0]] = 1
    else:
        dist = np.sqrt(X_xo[:, 0] ** 2 + X_xo[:, 1] ** 2)
        idx = np.argpartition(dist, 4)[:4]
        distp = np.prod(dist[idx])
        C[idx] = (distp / np.sum(distp / dist[idx])) / dist[idx]
    return C

vars={}
vars['A'] = AA2d
vars['B'] = BB2d
vars['Bv'] = BB2dv
vars['C_00_03'] = get_C((.02,.023))
vars['C_00__03'] = get_C((.02,.017))
vars['C_00_05'] = get_C((.02,.025))
vars['C_00__05'] = get_C((.02,.015))
vars['C_00_07'] = get_C((.02,.027))
vars['C_00__07'] = get_C((.02,.013))
vars['C_00_10'] = get_C((.02,.03))
vars['C_00__10'] = get_C((.02,.01))
vars['C_00_13'] = get_C((.02,.033))
vars['C_00__13'] = get_C((.02,.007))
vars['C_03_00'] = get_C((.023,.02))
vars['C__03_00'] = get_C((.017,.02))
vars['C_05_00'] = get_C((.025,.02))
vars['C__05_00'] = get_C((.015,.02))
vars['C_07_00'] = get_C((.027,.02))
vars['C__07_00'] = get_C((.013,.02))
vars['C_10_00'] = get_C((.03,.02))
vars['C__10_00'] = get_C((.01,.02))
vars['C_12_00'] = get_C((.032,.02))
vars['C__12_00'] = get_C((.008,.02))
vars['C_15_00'] = get_C((.035,.02))
vars['C_17_00'] = get_C((.037,.02))
vars['C_20_00'] = get_C((.04,.02))
vars['C_25_00'] = get_C((.045,.02))
vars['C_30_00'] = get_C((.05,.02))
vars['C_40_00'] = get_C((.06,.02))
import scipy.io
scipy.io.savemat('ss_vars_2d_gmsh_.mat', vars)