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
mesh = Mesh('../meshes/box_rev_symm.xml') # menos puntos para crear ss a usar en matlab
V = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

class DirichletBoundary2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > (0.12 - DOLFIN_EPS)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

# Pot=1500
ipx=.06
ipy=.0
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
p1=Point(0.06 + 0.014018,0,0)
p2=Point(0.06,0.002354,0)

# transitorio
u_prev = u1 #interpolate(u0, V)

u = TrialFunction(V)
v = TestFunction(V)

b = None
uf = Function(V)



t = 0.0
dt = .01
T = 150

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
ctrl1_d = control.sample_system(control.tf((2, 2), (1, 0)), dt)
ctrl1 = filter1.Filter1(ctrl1_d.num[0][0], ctrl1_d.den[0][0])
ctrl2_d = control.sample_system(control.tf((-5e-06, -1.5e-05), (3, 0)), dt)
ctrl2 = filter1.Filter1(ctrl2_d.num[0][0], ctrl2_d.den[0][0])

output_data = np.array([t,Pot,v_nom,u1(p1),u1(p2)])

ofile = File("T3Drsymm.pvd")
ofile << u_prev, t

exit(7)

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
p1=Point(0.06 - 0.014018,0.02,0)
# vdir = 1

while t <= T:
    t += dt
    Pt = Pot + ctrl1.step(er1)
    vt = (v_nom + ctrl2.step(er2)) * vdir
    a = rho*Cp*u*v*dx + dt*k*inner(nabla_grad(u), nabla_grad(v))*dx + dt*rho*Cp*inner(Constant((vt,0,0)), nabla_grad(u))*v*dx
    L = rho*Cp*u_prev*v*dx
    b = assemble(L, tensor=b)
    A = assemble(a)
    # bc.apply(A, b)
    delta = PointSource(V, inp, dt*Pt)
    delta.apply(b)
    solve(A, uf.vector(), b)
    u_prev.assign(uf)
    T1, T2 = u_prev(p1), u_prev(p2)
    output_data = np.vstack((output_data, [t, Pt, vt, T1, T2]))
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
        ofile << u_prev, t

ofile << u_prev, t
# print ("];", file=sys.stderr)
# np.save('output_data_vel_' + str(MPI.COMM_WORLD.Get_rank()), output_data)
# for op in o_points:
#     print("T(", *op, ",t=",T,"):", uf(*op), "K")
np.save('ctrl_data_test', output_data)
print("Tmax:", np.max(uf.vector().get_local()), "K")
print('T1:', T1, 'T2:', T2, 'Pt:', Pt, 'vt:', vt)
