
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import control
import weld1.filter1


df = pd.ExcelFile('C:\\tmp\\welding_log_2020-04-30_09-40-10.xlsx').parse('Weld Data')
tt = np.array(df['Time (s)'])
temp = np.array(df['Temperature (°C)'])

sim_wfs = np.zeros(tt.size)
sim_temp = np.zeros(tt.size)

Kp=7.5
Ki=10
Kd=0.5
Fd=.1
int_error = 0
# last_error = 0
d_state = 0
plant = control.tf(.3, (1/10, 1))
plant = control.tf(.3, (1/5, 1))
plant_d = control.sample_system(plant, .001)
fplant = weld1.filter1.Filter1((0,plant_d.num[0][0][0]),plant_d.den[0][0])

fig = plt.figure('Simulation animation')
ax = fig.add_subplot(111)
ax.clear()
ax.set_xlim((0,12.5))
ax.set_ylim((0,1200))
ax.plot(tt, temp-1100, label='temp')
lwfs, = ax.plot([], [], label='sim_wfs')
ltemp, = ax.plot([], [], label='sim_temp')
ax.legend()
ax.grid()

for i in range(1,tt.size):
    tsp = max(temp[i] - 1100, 0)
    # err = tsp - ftemp.step(sim_temp[i-1])
    error = tsp - sim_temp[i - 1]  # ftemp.step(sim_temp[i-1])
    # sim_wfs[i] = fctrl.step(err)
    int_error = int_error + (error * 0.001)
    d = Kd*Fd*error + d_state
    d_state = Kd*Fd*error + (1 - Fd * .001) * d
    # sim_wfs[i] = (Kp * error) + (Ki * int_error) + (Kd * (error - last_error) * 1000.0);
    sim_wfs[i] = (Kp * error) + (Ki * int_error) + d
    # last_error = error
    wfs1 = max(min(sim_wfs[i],1000), 200)
    if sim_wfs[i] > 1000:# or sim_wfs[i] > 200: # antiwindup
        int_error = int_error - (error * 0.001)
    sim_temp[i] = fplant.step(wfs1)
    if i % 100 == 0:
        lwfs.set_xdata(tt[:i])
        lwfs.set_ydata(sim_wfs[:i])
        ltemp.set_xdata(tt[:i])
        ltemp.set_ydata(sim_temp[:i])
        fig.canvas.draw()
        fig.canvas.flush_events()

# plt.figure('Plc Simulation')
# plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, sim_wfs, label='sim_wfs')
# plt.plot(tt, sim_temp, label='sim_temp')
# plt.legend()
# plt.grid()
# plt.ylim((0,1200))