#  docker run -ti -v C:/src/sympy:/home/fenics/shared fenics1
from dolfin import *
import numpy as np


mesh = Mesh('../meshes/box_reverse.xml') # menos puntos para crear ss a usar en matlab
V = FunctionSpace(mesh, "CG", 1)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# Define boundary condition
u0 = Constant(0.0)
# bc = DirichletBC(V, u0, DirichletBoundary())
bc = DirichletBC(V, u0, 'on_boundary && near(x[0], 0, 1E-14)')
bc2 = DirichletBC(V, u0, 'on_boundary && near(x[0], 0.12, 1E-14)')

ipx=.06
ipy=.02
inp = Point(ipx,ipy,.0)

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

p1=Point(0.06 + 0.014534,0.02,0)
p2=Point(0.06,0.022353,0)

class matrixes:
    AAA = None
    BB = None
    BBv = None
    TT0 = None
    TT0r = None


mtxs = matrixes()

save_ss = True

def plant(k_=60,cp_=2719,pot_=7500,v_=.005):

    # system from references isotherms
    k = Constant(k_)
    rho = 7870
    Cp = Constant(cp_)
    Pot = pot_
    v_nom = v_
    vel = Constant((v_nom,0,0))

    a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
    L1 = Constant(0)*v*dx

    A1, b1 = assemble_system(a1, L1, bc)
    delta = PointSource(V, inp, Pot)
    delta.apply(b1)

    u1 = Function(V)

    # Stationary
    solve(A1, u1.vector(), b1)
    # print ("max: ", max(u1.vector()),"(K), T(x=0.2): ", u1.vector()[500]) # 137.15

    print("k,cp,pot,v:",k_,cp_,pot_,v_)
    print("T1:",str(u1(p1)),"T2:",str(u1(p2)))

    # save_ss = False
    if save_ss:
        K, f = assemble_system(a1, L1, bc)
        delta.apply(f)
        a = rho*Cp*u*v*dx
        C = assemble(a)
        # bc.apply(C)

        # gT0 = project(T0.dx(0), V)

        a = rho*Cp*inner(Constant((1,0,0)), nabla_grad(u))*v*dx
        Kv = assemble(a)
        # bc.apply(Kv)

        # print('About to solve AA')
        AA3d = - np.linalg.solve(C.array(), K.array())
        BB3d = np.array([np.linalg.solve(C.array(), f.get_local() / Pot)]).T
        BB3dv = - np.array([np.linalg.solve(C.array(), Kv.array() @ u1.vector().get_local())]).T

        vel = Constant((-v_, 0, 0))
        a_1 = k * inner(nabla_grad(u), nabla_grad(v)) * dx + rho * Cp * inner(vel, nabla_grad(u)) * v * dx
        A_1, b_1 = assemble_system(a_1, L1, bc2)
        delta_1 = PointSource(V, inp, Pot)
        delta_1.apply(b_1)
        u_1 = Function(V)

        # Stationary
        solve(A_1, u_1.vector(), b_1)

        if mtxs.AAA is None:
            mtxs.AAA = AA3d
            mtxs.BB = BB3d
            mtxs.BBv = BB3dv
            mtxs.TT0 = np.array([u1.vector().get_local()]).T
            mtxs.TT0r = np.array([u_1.vector().get_local()]).T
        else:
            mtxs.AAA = np.dstack((mtxs.AAA, AA3d))
            mtxs.BB = np.dstack((mtxs.BB, BB3d))
            mtxs.BBv = np.dstack((mtxs.BBv, BB3dv))
            mtxs.TT0 = np.dstack((mtxs.TT0, np.array([u1.vector().get_local()]).T))
            mtxs.TT0r = np.dstack((mtxs.TT0r, np.array([u_1.vector().get_local()]).T))

        # print('About to store XX, AA, ...')
        # np.save('XX3d', dof_coordinates)
        # np.save('AA3d', AA3d)
        # np.save('BB3d', BB3d)
        # np.save('BB3dv', BB3dv)
        # np.save('T03d', np.array([u1.vector().get_local()]).T)
        # np.save('p1', np.array((p1[0],p1[1],p1[2])))
        # np.save('p2', np.array((p2[0], p2[1], p2[2])))
        # print('mtxs.AAA.shape', mtxs.AAA.shape)


for k in (50,60,70):
    for cp in (2500,2719,3000):
        for vel in (.0045, .005, .0055):
            plant(k_=k, cp_=cp, v_=vel)


if save_ss:
    np.save('AAA', mtxs.AAA)
    np.save('BB', mtxs.BB)
    np.save('BBv', mtxs.BBv)
    np.save('TT0', mtxs.TT0)
    np.save('TT0r', mtxs.TT0r)