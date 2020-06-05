from dolfin import *
import numpy as np

# Create mesh and function space
mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 64, 32, 8)
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

k = Constant(24)
rho = 7925
Cp = 460

vel = Constant((0.002,0,0))

a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

Pot=1000
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


from scipy.sparse import csr_matrix, save_npz
from scipy.sparse.linalg import spsolve
C_mat = as_backend_type(C).mat()
K_mat = as_backend_type(K).mat()
Kv_mat = as_backend_type(Kv).mat()
C_sparray = csr_matrix(C_mat.getValuesCSR()[::-1], shape = C_mat.size)
K_sparray = csr_matrix(K_mat.getValuesCSR()[::-1], shape = K_mat.size)
Kv_sparray = csr_matrix(Kv_mat.getValuesCSR()[::-1], shape = Kv_mat.size)

save_npz('C_sp', C_sparray)
save_npz('K_sp', K_sparray)

exit(0)

from scipy.sparse import load_npz, save_npz
from scipy.sparse.linalg import spsolve
C_sparray = load_npz('C_sp.npz')
K_sparray = load_npz('K_sp.npz')
print('Converting sparse csr to csc')
C_sp = C_sparray.tocsc()
K_sp = K_sparray.tocsc()
print('Solving C\\K')
AA3d = -spsolve(C_sp, K_sp)
print('Saving')
save_npz('AA3d_sp', AA3d)
BB3d = np.array([spsolve(C_sparray, f.get_local() / Pot)]).T
BB3dv = -np.array([spsolve(C_sparray, K_sparray @ T0.vector().get_local())]).T

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

np.save('XX3d', dof_coordinates)
# np.save('AA3d', AA3d)
save_npz('AA3d', AA3d)
np.save('BB3d', BB3d)
np.save('BB3dv', BB3dv)