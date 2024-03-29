#  docker run -ti -v C:/src/sympy:/home/fenics/shared fenics1

from dolfin import *
import numpy as np
import control
import filter1


# Create mesh and function space
# mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 320, 160, 40)
# mesh = Mesh('meshes/box.xml') # mas preciso
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

# bcf = Final()
# for x in mesh.coordinates():
#     if bcf.inside(x, True): print('%s is on x = 0' % x)

# exit(0)
# Pot=1500
ipx=.06
ipy=.02
inp = Point(ipx,ipy,.0)

v_nom = 0.005

# system from references isotherms
k = Constant(60)
rho = 7870
Cp = 2719
Pot = 7500
v_nom = 0.005
vel = Constant((v_nom,0,0))


# a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*vel*u.dx(0)*v*dx + P*h*inner(u,v)*dx
# L1 = Constant(0)*v*dx + P*h*Tinf*v*dx
a1 = k*inner(nabla_grad(u), nabla_grad(v))*dx + rho*Cp*inner(vel, nabla_grad(u))*v*dx
L1 = Constant(0)*v*dx

A1, b1 = assemble_system(a1, L1, bc)
delta = PointSource(V, inp, Pot)
delta.apply(b1)

u1 = Function(V)

# Stationary
solve(A1, u1.vector(), b1)
# print ("max: ", max(u1.vector()),"(K), T(x=0.2): ", u1.vector()[500]) # 137.15

n = V.dim()
d = mesh.geometry().dim()
dof_coordinates = V.tabulate_dof_coordinates()
dof_coordinates.resize((n, d))

p1=Point(0.06 + 0.014534,0.02,0)
p2=Point(0.06,0.022353,0)

T1ref=1200
T2ref=1200
# # box
# p1=Point(0.036639,0.02,0)
# p2=Point(0.02,0.022333,0)

# box2
p1=Point(0.06 + 0.014534,0.02,0)
p2=Point(0.06,0.022353,0)

# transitorio
u_prev = u1 #interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

b = None
uf = Function(V)



t = 0.0
dt = .01
T = 500

# Controllers Jorge different QFT refs 1
ctrl1_d = control.sample_system(control.tf((1.8,), (1, 0)), dt) #Gc_pot = tf(1.8,[1 0]);
ctrl1 = filter1.Filter1((0, ctrl1_d.num[0][0][0]), ctrl1_d.den[0][0])
ctrl2_d = control.sample_system(control.tf((-9.3e-6,), (1, 0)), dt) #Gc_vel = tf(-9.3e-6,[1 0]);
ctrl2 = filter1.Filter1((0, ctrl2_d.num[0][0][0]), ctrl2_d.den[0][0])

output_data = np.array([t,Pot,v_nom,u1(p1),u1(p2)])

save_pvd = False
if save_pvd:
    ofile = File("T3Dr_c2.pvd")
    ofile << u_prev, t

tn=1
kt=10
T1,T2=u_prev(p1),u_prev(p2)
er1,er2=T1ref - T1,T2ref - T2
print("T1:", T1, "T2:", T2)
stable=0

# Perturbations
# k = Constant(20)
# # Step
T1ref += 1000
# T2ref += -100
# upert = u_prev.vector().get_local() + 100
# u_prev.vector().set_local(upert)
# vdir = -1.0
# # p1=Point(0.005982,0.02,0)
# p1=Point(0.06 - 0.014534,0.02,0)
vdir = 1

while t <= T:
    t += dt
    # Pt = Pot + ctrl1.step(er1)
    # vt = (v_nom + ctrl2.step(er2)) * vdir
    Pt = Pot
    vt = v_nom * 1.1
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*inner(Constant((vt,0,0)), nabla_grad(u))*v*dx
    L = rho*Cp*u_prev*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    bc.apply(A, b)
    # bc2.apply(A, b)
    delta = PointSource(V, inp, dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    u_prev.assign(uf)
    T1, T2 = u_prev(p1), u_prev(p2)
    output_data = np.vstack((output_data, [t, Pt, vt, T1, T2]))

    er1, er2 = T1ref - T1, T2ref - T2
    if (er1**2+er2**2) < .01:
        stable += 1
        if stable > 20:
            print("Stable output t=", round(t, 3), 'of', T, 'seconds T1:', T1, 'T2:', T2)
            break
    else:
        stable = 0
    if t*kt >= tn:
        print("t=", round(t,3), 'of', T, 'seconds T1:', T1, 'T2:', T2)
        # print("er1:", er1, "er2:", er2, "norm(er):", (er1**2+er2**2))
        tn += 1
        if t > 5 and kt == 10:
            kt = 1
            tn = int(t*kt) + 1
        if t > 25 and kt == 1:
            kt = .1
            tn = int(t*kt) + 1
        if save_pvd:
            ofile << u_prev, t

if save_pvd:
    ofile << u_prev, t
# print ("];", file=sys.stderr)
# np.save('output_data_vel_' + str(MPI.COMM_WORLD.Get_rank()), output_data)
# for op in o_points:
#     print("T(", *op, ",t=",T,"):", uf(*op), "K")
np.save('ctrl_steps_data2', output_data)
print("Tmax:", np.max(uf.vector().get_local()), "K")
print('T1:', T1, 'T2:', T2, 'Pt:', Pt, 'vt:', vt)
