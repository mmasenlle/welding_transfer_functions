
import matplotlib.pyplot as plt
import numpy as np
import control

import fem3d_2ss, weld1.filter1
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/16x8x2/')
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/1500/')
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/40x20x5/')
# fem2ss_2 = fem3d_2ss.Fem3d_fenics('data/v5/32x16x4/1500/k30/')
# fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/irregular/')
fem2ss = fem3d_2ss.Fem3d_fenics('data/v5/irregular2/')
fem2ss_2 = fem3d_2ss.Fem3d_fenics('data/v5/irregular2/k30/')
# fem2ss.A.shape


inp = np.array((.02,.02,0))
p1 = (.01,0,0)
p2 = (0,0.005,0)
plantp_d = control.sample_system(fem2ss.get_ss(inp + p1), .001)
plantv_d = control.sample_system(fem2ss.get_ss_v(inp + p2), .001)
# plant_d = plantp_d

T = fem2ss.Teq

control_on = True
perturb_on = False

# Controllers
ctrl1_d = control.sample_system(control.tf(1, (1, 0)), .001)
ctrl1 = weld1.filter1.Filter1((0,ctrl1_d.num[0][0][0]),ctrl1_d.den[0][0])
ctrl2_d = control.sample_system(control.tf(-0.0002, (1, 0)), .001)
ctrl2 = weld1.filter1.Filter1((0,ctrl2_d.num[0][0][0]),ctrl2_d.den[0][0])

# Set points
power = 1500.0
delta_speed = 0.0
T1ref = 1000
T2ref = 400
power0 = power

# Simulation. T1, T2
fig = plt.figure('Temperatures')
dt=.001
tt=np.arange(0, 10 + dt, dt)
T1 = np.zeros(tt.size)
T2 = np.zeros(tt.size)
ctrl_power = np.zeros(tt.size)
ctrl_speed = np.zeros(tt.size)
labels_unit_scale=(('T1','(ºC)',2000),('Power','(w)',3000),('T2','(ºC)',1000),('Speed','(m/s)',0.01))
ydata = (T1, ctrl_power, T2, ctrl_speed)
axs,line,n=[],[],len(ydata)
for i in range(n):
    ax = fig.add_subplot(100 * n + 11 + i)
    lin, = ax.plot(tt, ydata[i], label=labels_unit_scale[i][0])
    ax.set_ylabel(labels_unit_scale[i][1])
    ax.set_ylim((0, labels_unit_scale[i][2]))
    ax.legend()
    ax.grid()
    axs.append(ax)
    line.append(lin)
axs[0].set_title('Control temperature (0 s)')
axs[3].set_xlabel('(s)')

for i in range(tt.size):
    T = plantp_d.A @ T + plantp_d.B * power + plantv_d.B * delta_speed
    T1[i] = plantp_d.C @ T
    T2[i] = plantv_d.C @ T
    if control_on and i > 0:
        power = power0 + ctrl1.step(T1ref - T1[i])
        delta_speed = ctrl2.step(T2ref - T2[i])
    ctrl_power[i] = power
    ctrl_speed[i] = delta_speed + .005
    if perturb_on and i > 5000:
        perturb_on = False
        plantp_d = control.sample_system(fem2ss_2.get_ss(inp + p1), .001)
        plantv_d = control.sample_system(fem2ss_2.get_ss_v(inp + p2), .001)
    if i % 100 == 0:
        for k in range(n):
            line[k].set_ydata(ydata[k])
        axs[0].set_title('T (t='+str(round(i/1000, 2))+' s)')
        fig.canvas.draw()
        fig.canvas.flush_events()



# Simulation. Temp profiles along x and y
nds1 = np.argwhere((fem2ss.X[:,1]==0.02)&(fem2ss.X[:,2]==0))
nds2 = np.argwhere((fem2ss.X[:,0]==0.02)&(fem2ss.X[:,2]==0))
# nds1 = np.argwhere((np.abs(fem2ss.X[:,1]-0.02) < .001)&(fem2ss.X[:,2]==0))
# nds2 = np.argwhere((np.abs(fem2ss.X[:,0]-0.02) < .001)&(fem2ss.X[:,2]==0))

fig = plt.figure('Temp profiles 2')
ax1 = fig.add_subplot(211)
ax1.clear()
ax1.set_xlim((0,.08))
ax1.set_ylim((0,3000))
line1, = ax1.plot(fem2ss.X[nds1,0], T[nds1,0], '.', label='T along x')
line2, = ax1.plot(fem2ss.X[nds2,1], T[nds2,0], '.', label='T along y')
ax1.plot(inp[0] + p1[0], T1ref, marker='+', color='C0')
ax1.plot(inp[1] + p2[1], T2ref, marker='+', color='C1')
ax1.set_title('T (t=0 s)')
ax1.legend()
ax1.grid()
ax1.set_xlabel('(m)')
ax1.set_ylabel('(C)')
dt=.001
tt=np.arange(0, 10 + dt, dt)


ax2 = fig.add_subplot(212)
ax2.set_ylim((0,3000))
ax2.set_xlabel('(s)')
ax2.set_ylabel('(w)')
ctrl_power = np.zeros(tt.size)
ctrl_speed = np.zeros(tt.size)
line3, = ax2.plot(tt, ctrl_power, label='power')
line4, = ax2.plot(tt, ctrl_speed, label='speed')
ax2.legend()
power0 = power
for i in range(tt.size):
    T = plantp_d.A @ T + plantp_d.B * power + plantv_d.B * delta_speed
    if control_on and i > 0:
        power = power0 + ctrl1.step(T1ref - (plantp_d.C @ T))
        delta_speed = ctrl2.step((T2ref - plantv_d.C @ T))
    ctrl_power[i] = power
    ctrl_speed[i] = (delta_speed + .005) * 100000
    if perturb_on and i > 5000:
        perturb_on = False
        plantp_d = control.sample_system(fem2ss_2.get_ss(inp + p1), .001)
        plantv_d = control.sample_system(fem2ss_2.get_ss_v(inp + p2), .001)
    if i % 100 == 0:
        line1.set_ydata(T[nds1,0])
        line2.set_ydata(T[nds2, 0])
        ax1.set_title('T (t='+str(round(i/1000, 2))+' s)')
        line3.set_ydata(ctrl_power)
        line4.set_ydata(ctrl_speed)
        fig.canvas.draw()
        fig.canvas.flush_events()
print(power,delta_speed)


# Surface
# power = 500
nds4 = np.argwhere((fem2ss.X[:,2]==0))
# T = np.zeros(fem2ss.B.shape)
from scipy.interpolate import interp2d

x_coords = np.arange(min(fem2ss.X[nds4,0]), max(fem2ss.X[nds4,0]), .0001)
y_coords = np.arange(min(fem2ss.X[nds4,1]), max(fem2ss.X[nds4,1]), .0001)

f = interp2d(fem2ss.X[nds4,0], fem2ss.X[nds4,1], T[nds4,0], kind="linear")
Z = f(x_coords,y_coords)

fig2 = plt.figure('Temp surface')
img = plt.imshow(Z,
           extent=[min(x_coords),max(x_coords),min(y_coords),max(y_coords)],
           origin="lower", vmin=0, vmax=3000)
plt.title('T of surface')
plt.colorbar()
plt.plot(inp[0] + p1[0], inp[1] + p1[1], marker='+', color='k')
plt.plot(inp[0] + p2[0], inp[1] + p2[1], marker='+', color='k')
plt.ylabel('m')

plt.figure('Profiles temp')
plt.plot(x_coords, Z[200,:], label='T along x')
plt.plot(y_coords, Z[:,200], label='T along y')


# for i in range(20001):
#     T = plantp_d.A @ T + plantp_d.B * power + plantv_d.B * delta_speed
#     if i % 500 == 0:
#         plt.title('T (t='+str(round(i/1000, 2))+' s)')
#         f = interp2d(fem2ss.X[nds2, 0], fem2ss.X[nds2, 1], T[nds2, 0], kind="linear")
#         Z = f(x_coords, y_coords)
#         img.set_data(Z)
#         # img = plt.imshow(Z,
#         #              extent=[min(x_coords), max(x_coords), min(y_coords), max(y_coords)],
#         #              origin="lower", vmin = 0, vmax = 3000)
#         fig2.canvas.draw()
#         fig2.canvas.flush_events()



# Section
nds3 = np.argwhere((fem2ss.X[:,1]==0.02))
# T = np.zeros(fem2ss.B.shape)
from scipy.interpolate import interp2d

x_coords = np.arange(min(fem2ss.X[nds3,0]), max(fem2ss.X[nds3,0]), .0001)
y_coords = np.arange(min(fem2ss.X[nds3,2]), max(fem2ss.X[nds3,2]), .0001)

f = interp2d(fem2ss.X[nds3,0], fem2ss.X[nds3,2], T[nds3,0], kind="linear")
Z = f(x_coords,y_coords)

fig3 = plt.figure('Temp section')
img = plt.imshow(Z,
           extent=[min(x_coords),max(x_coords),min(y_coords),max(y_coords)],
           origin="upper", vmin=0, vmax=3000)
plt.title('T of section')
plt.colorbar()
plt.xlabel('m')
plt.ylabel('m')

# for i in range(20001):
#     T = plantp_d.A @ T + plantp_d.B * power + plantv_d.B * delta_speed
#     if i % 500 == 0:
#         plt.title('T (t='+str(round(i/1000, 2))+' s)')
#         f = interp2d(fem2ss.X[nds3, 0], fem2ss.X[nds3, 2], T[nds3, 0], kind="linear")
#         Z = f(x_coords, y_coords)
#         img.set_data(Z)
#         # img = plt.imshow(Z,
#         #              extent=[min(x_coords), max(x_coords), min(y_coords), max(y_coords)],
#         #              origin="lower", vmin = 0, vmax = 3000)
#         fig2.canvas.draw()
#         fig2.canvas.flush_events()

