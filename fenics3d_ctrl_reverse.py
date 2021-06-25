#  docker run -ti -v C:/src/sympy:/home/fenics/shared fenics1
import sys
from dolfin import *
import numpy as np
import control
import filter1
from mpi4py import MPI

def on_message(mqttc, obj, msg):
    print(msg.topic + " " + str(msg.qos) + " " + str(msg.payload))

mqtt_client=None
b_mqtt = True
if b_mqtt:
    try:
        import paho.mqtt.client as mqtt
        mqtt_client = mqtt.Client()
        mqtt_client.on_message = on_message
        mqtt_client.connect('172.24.1.45')
        mqtt_client.subscribe("ctrl3d_cmd", 0)
        mqtt_client.loop_start()
    except:
        print("Cant connect mqtt broker")
        b_mqtt = False

# Create mesh and function space
# mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(.08, .04, .01), 320, 160, 40)
# mesh = Mesh('meshes/box.xml') # mas preciso
mesh = Mesh('meshes/box_reverse.xml') # menos puntos para crear ss a usar en matlab
V = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
# class DirichletBoundary(SubDomain):
#     def inside(self, x, on_boundary):
#         return on_boundary and x[0] < DOLFIN_EPS
#
# class Final(SubDomain):
#     tol = 1E-14
#     def inside(self, x, on_boundary):
#         return on_boundary and near(x[0], 0.12, self.tol)

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

# # Equilibrium
# Pot=1481.6279283811557
# v_nom=0.004253872403429254
#
#
# k = Constant(24)
# rho = Constant(7925)
# #Pot=2
# #k = Constant(0.000024)
# #rho = Constant(0.007925)
# Cp = Constant(460)


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

save_ss = True
save_ss = False
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

    print('About to solve AA')
    AA3d = - np.linalg.solve(C.array(), K.array())
    BB3d = np.array([np.linalg.solve(C.array(), f.get_local() / Pot)]).T
    BB3dv = - np.array([np.linalg.solve(C.array(), Kv.array() @ u1.vector().get_local())]).T

    print('About to store XX, AA, ...')
    np.save('XX3d', dof_coordinates)
    np.save('AA3d', AA3d)
    np.save('BB3d', BB3d)
    np.save('BB3dv', BB3dv)
    np.save('T03d', np.array([u1.vector().get_local()]).T)
    np.save('p1', np.array((p1[0],p1[1],p1[2])))
    np.save('p2', np.array((p2[0], p2[1], p2[2])))

# T1ref=1000
# T2ref=400
# p1=Point(0.03,0.02,0)
# p2=Point(0.02,0.025,0)

# isotherms 1200
# T(0.034018,0.02,0): 1199.997701
# T(0.02,0.022354,0): 1198.289010
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
T = 50

# # Controllers
# ctrl1_d = control.sample_system(control.tf(2, (1, 0)), dt)
# ctrl1 = filter1.Filter1((0,ctrl1_d.num[0][0][0]),ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf(-0.0001, (1, 0)), dt)
# ctrl2 = filter1.Filter1((0,ctrl2_d.num[0][0][0]),ctrl2_d.den[0][0])
#
# # Controllers Jorge MF=30
# ctrl1_d = control.sample_system(control.tf((2, 2.6), (1.3, 0)), dt)
# ctrl1 = filter1.Filter1(ctrl1_d.num[0][0],ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf((-0.0001021,-0.0004546), (4.453, 45.46, 0)), dt)
# ctrl2 = filter1.Filter1((0,ctrl2_d.num[0][0][0],ctrl2_d.num[0][0][1]),ctrl2_d.den[0][0])
#
# # Controllers Jorge MF=40
# ctrl1_d = control.sample_system(control.tf((1, 3), (3, 0)), dt)
# ctrl1 = filter1.Filter1(ctrl1_d.num[0][0],ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf(-5e-06, (1, 1, 0)), dt)
# ctrl2 = filter1.Filter1((0,ctrl2_d.num[0][0][0],ctrl2_d.num[0][0][1]),ctrl2_d.den[0][0])

# Controllers Jorge different bandwithes
# ctrl1_d = control.sample_system(control.tf((2, 2), (1, 0)), dt)
# ctrl1_d = control.sample_system(control.tf((1, 1), (1, 0)), dt)
# ctrl1 = filter1.Filter1(ctrl1_d.num[0][0], ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf((-5e-06, -1.5e-05), (3, 0)), dt)
# ctrl2 = filter1.Filter1(ctrl2_d.num[0][0], ctrl2_d.den[0][0])

# Controllers Jorge different QFT refs 1
# ctrl1_d = control.sample_system(control.tf((1.8,), (1, 0)), dt) #Gc_pot = tf(1.8,[1 0]);
# ctrl1 = filter1.Filter1((0, ctrl1_d.num[0][0][0]), ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf((-9.3e-6,), (1, 0)), dt) #Gc_vel = tf(-9.3e-6,[1 0]);
# ctrl2 = filter1.Filter1((0, ctrl2_d.num[0][0][0]), ctrl2_d.den[0][0])

# Controllers test linealidad QFT refs
# ctrl1_d = control.sample_system(control.tf((.2,), (1, 0)), dt)
# ctrl1 = filter1.Filter1((0, ctrl1_d.num[0][0][0]), ctrl1_d.den[0][0])
# ctrl2_d = control.sample_system(control.tf((-1e-6,), (1, 0)), dt)
# ctrl2 = filter1.Filter1((0, ctrl2_d.num[0][0][0]), ctrl2_d.den[0][0])

# Controllers and filters QFT manu 2
# F = [-50000/(s/0.1+1)/(s/1+1); %tf(0);
#     0.005*(s/2+1)/(s/1+1)/(s/4+1)];
# G = [1/s^1*(s/1+1)*(s/2+1)/(s/10+1) tf(0);
#     tf(0) -4e-05/s^1*(s/5+1)/(s/40+1)];
# ctrl1_d = control.sample_system(control.tf((10,30,20), (2, 20, 0)), dt)
# ctrl1 = filter1.Filter1(ctrl1_d.num[0][0], ctrl1_d.den[0][0])
# # ctrl2_d = control.sample_system(control.tf((-0.0016,-0.008), (5, 200, 0)), dt)
# ctrl2_d = control.sample_system(control.tf((-7.199e-05,-0.0003924,-0.0001622), (2.253, 92.13, 80.14, 0)), dt)
# ctrl2 = filter1.Filter1(ctrl2_d.num[0][0], ctrl2_d.den[0][0])
# # flt1_d = control.sample_system(control.tf((-5000,), (1, 1.1, 0.1)), dt)
# flt1_d = control.sample_system(control.tf((-60000,), (1, 10.1, 1)), dt)
# flt1 = filter1.Filter1(flt1_d.num[0][0], flt1_d.den[0][0])
# # flt2_d = control.sample_system(control.tf((.02, .04), (2, 10, 8)), dt)
# flt2_d = control.sample_system(control.tf((.02,), (1, 5, 4)), dt)
# flt2 = filter1.Filter1(flt2_d.num[0][0], flt2_d.den[0][0])

# Controllers and filters QFT manu 3
ctrl1_d = control.sample_system(control.tf((10,), (1, 10, 0)), dt)
ctrl1 = filter1.Filter1(ctrl1_d.num[0][0], ctrl1_d.den[0][0])
ctrl2_d = control.sample_system(control.tf((-0.0002,-0.00065,-0.000375), (1.875, 189.4, 187.5, 0)), dt)
ctrl2 = filter1.Filter1(ctrl2_d.num[0][0], ctrl2_d.den[0][0])
flt1_d = control.sample_system(control.tf((-50000,), (1, 10.1, 1)), dt)
flt1 = filter1.Filter1(flt1_d.num[0][0], flt1_d.den[0][0])
flt2_d = control.sample_system(control.tf((.0114,), (1, 3.8)), dt)
flt2 = filter1.Filter1(flt2_d.num[0][0], flt2_d.den[0][0])


output_data = np.array([t,Pot,v_nom,u1(p1),u1(p2),0,0])

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
# T1ref += -100
# T2ref += -100
# upert = u_prev.vector().get_local() + 100
# u_prev.vector().set_local(upert)
vdir = -1.0
# p1=Point(0.005982,0.02,0)
p1=Point(0.06 - 0.014534,0.02,0)
# vdir = 1
impul = 1/dt

while t <= T:
    t += dt
    Pff = flt1.step(impul)
    vff = flt2.step(impul)
    Pt = Pot + ctrl1.step(er1) + Pff
    vt = (v_nom + ctrl2.step(er2) + vff) * vdir
    impul = 0
    # Pt = Pot
    # vt = v_nom * vdir
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*inner(Constant((vt,0,0)), nabla_grad(u))*v*dx
    L = rho*Cp*u_prev*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    # bc.apply(A, b)
    bc2.apply(A, b)
    delta = PointSource(V, inp, dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    u_prev.assign(uf)
    T1, T2 = u_prev(p1), u_prev(p2)
    output_data = np.vstack((output_data, [t, Pt, vt, T1, T2, Pff, vff]))
    if b_mqtt:
        mqtt_client.publish('ctrl3d', output_data[-1,:].tobytes())
    er1, er2 = T1ref - T1, T2ref - T2
    if (er1**2+er2**2) < .01:
        stable += 1
        if stable > 20:
            print("Stable output t=", round(t, 3), 'of', T, 'seconds T1:', T1, 'T2:', T2)
            break
    else:
        stable = 0
    if t*kt >= tn:
        print("t=", round(t,3), 'of', T, 'seconds')
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
np.save('ctrl_qft_manu3', output_data)
print("Tmax:", np.max(uf.vector().get_local()), "K")
print('T1:', T1, 'T2:', T2, 'Pt:', Pt, 'vt:', vt)

if save_ss:
    np.save('T03dr', np.array([uf.vector().get_local()]).T)
    # vel = Constant((-v_nom, 0, 0))
    # a1 = k * inner(nabla_grad(u), nabla_grad(v)) * dx + rho * Cp * inner(vel, nabla_grad(u)) * v * dx
    # L1 = Constant(0) * v * dx
    #
    # K, f = assemble_system(a1, L1, bc2)
    # delta.apply(f)
    # T0 = Function(V)
    # solve(K, T0.vector(), f)
    #
    # a = rho * Cp * u * v * dx
    # C = assemble(a)
    # # bc.apply(C)
    #
    # # gT0 = project(T0.dx(0), V)
    #
    # a = rho * Cp * inner(Constant((1, 0, 0)), nabla_grad(u)) * v * dx
    # Kv = assemble(a)
    # # bc.apply(Kv)
    #
    # print('About to solve AA')
    # AA3d = - np.linalg.solve(C.array(), K.array())
    # BB3d = np.array([np.linalg.solve(C.array(), f.get_local() / Pot)]).T
    # BB3dv = - np.array([np.linalg.solve(C.array(), Kv.array() @ T0.vector().get_local())]).T
    #
    # print('About to store XX, AA, ...')
    # np.save('XX3dr', dof_coordinates)
    # np.save('AA3dr', AA3d)
    # np.save('BB3dr', BB3d)
    # np.save('BB3dvr', BB3dv)
    # np.save('T03dr', np.array([T0.vector().get_local()]).T)
    # np.save('p1r', np.array((p1[0], p1[1], p1[2])))
    # np.save('p2r', np.array((p2[0], p2[1], p2[2])))
