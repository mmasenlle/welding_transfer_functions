
import matplotlib.pyplot as plt
import numpy as np
import control

import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics()

plant = fem2ss.get_ss((.02,.025,0))
plant_d = control.sample_system(plant, .001)

nds1 = np.argwhere((fem2ss.X[:,1]==0.02)&(fem2ss.X[:,2]==0))
T = np.zeros(fem2ss.B.shape)
# plt.plot(fem2ss.X[nds1,0], T[nds1,0])

fig = plt.figure('Temp profile')
ax = fig.add_subplot(111)
ax.clear()
ax.set_xlim((0,.08))
ax.set_ylim((0,3000))
line, = ax.plot(fem2ss.X[nds1,0], T[nds1,0], label='temp')
ax.set_title('T (t=0 s)')
ax.legend()
ax.grid()

power = 500

for i in range(10001):
    T = plant_d.A @ T + plant_d.B * power
    if i % 50 == 0:
        line.set_ydata(T[nds1,0])
        ax.set_title('T (t='+str(round(i/1000, 2))+' s)')
        fig.canvas.draw()
        fig.canvas.flush_events()

# Surface
nds2 = np.argwhere((fem2ss.X[:,2]==0))
T = np.zeros(fem2ss.B.shape)
from scipy.interpolate import interp2d

x_coords = np.arange(min(fem2ss.X[nds2,0]), max(fem2ss.X[nds2,0]), .001)
y_coords = np.arange(min(fem2ss.X[nds2,1]), max(fem2ss.X[nds2,1]), .001)

f = interp2d(fem2ss.X[nds2,0], fem2ss.X[nds2,1], T[nds2,0], kind="linear")
Z = f(x_coords,y_coords)

fig2 = plt.figure('Temp surface')
img = plt.imshow(Z,
           extent=[min(x_coords),max(x_coords),min(y_coords),max(y_coords)],
           origin="lower", vmin = 0, vmax = 3000)
plt.title('T (t=0 s)')
plt.colorbar()

for i in range(20001):
    T = plant_d.A @ T + plant_d.B * power
    if i % 100 == 0:
        plt.title('T (t='+str(round(i/1000, 2))+' s)')
        f = interp2d(fem2ss.X[nds2, 0], fem2ss.X[nds2, 1], T[nds2, 0], kind="linear")
        Z = f(x_coords, y_coords)
        #     img.set_data(Z)
        img = plt.imshow(Z,
                     extent=[min(x_coords), max(x_coords), min(y_coords), max(y_coords)],
                     origin="lower", vmin = 0, vmax = 3000)
        fig2.canvas.draw()
        fig2.canvas.flush_events()


