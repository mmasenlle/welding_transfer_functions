import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


df = pd.ExcelFile('C:\\tmp\\welding_log_2020-04-30_09-40-10.xlsx').parse('Weld Data')

tt = np.array(df['Time (s)'])
dt = np.diff(tt)
max(dt)
min(dt)
temp = np.array(df['Temperature (°C)'])
wfs = np.array(df['Wire Speed (m/min)'])
swfs = np.array(df['Wire Speed S (m/min)'])


plt.figure('Wire Speed')
plt.plot(tt,wfs*100, label='fronius')
plt.plot(tt,swfs*100, label='sensor')
plt.plot(tt,temp, label='temp')
plt.legend()

vv = np.array(df['Voltage (V)'])
aa = np.array(df['Current (A)'])
saa = np.array(df['Current S (A)'])

plt.figure('Power')
# plt.plot(tt,vv, label='voltage')
# plt.plot(tt,aa, label='current')
# plt.plot(tt,saa, label='current s')
# plt.plot(tt,aa * vv, label='power (w)')
plt.plot(tt,saa * vv, label='power (w)')
plt.plot(tt, 2*temp, label='2 x Temp (C)')
plt.plot(tt,swfs*100, label='wfs')
plt.legend()


dt = tt[1]-tt[0]
N = len(tt)
freqs = np.arange(N)/(N*dt)
N2=int(N/2)

fft_t = np.fft.fft(temp)
plt.figure('fft')
plt.plot(freqs[:N2], np.abs(fft_t[:N2]))
plt.ylim([0,10000])

from scipy import signal
b, a = signal.butter(3, 50/500, 'low')
temp_filtered_50 = signal.filtfilt(b, a, temp)
plt.figure('Temp')
plt.plot(tt,temp, label='data')
plt.plot(tt, temp_filtered_50, label='temp_filtered_50')
b, a = signal.butter(3, 10/500, 'low')
temp_filtered_10 = signal.filtfilt(b, a, temp)
plt.plot(tt, temp_filtered_10, label='temp_filtered_10')
b, a = signal.butter(3, 2/500, 'low')
temp_filtered_2 = signal.filtfilt(b, a, temp)
plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
plt.legend()

# import control
# mig = control.tf(1, (1/10 + 1))
plant = signal.lti(.3, (1/10, 1))
t,y = signal.step(plant, T=tt)
plant2 = signal.lti(.3, (1/1, 1))
t2,y2 = signal.step(plant2, T=tt)
plt.figure('Plant')
plt.plot(t, 700*y, label='plant')
plt.plot(t2, 700*y2, label='plant2')
plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
plt.ylim((0,400))
plt.legend()


import control
import control.matlab as mt
Kp=10
Ki=60
Kd=0.2
ctrl = control.tf(Kp, 1) + control.tf(Ki, (1,0)) + control.tf((Kd,0), 1)
# control.pzmap(ctrl)
plant = control.tf(.3, (1/10, 1))
flt = control.tf(b, a, .001)
# t1,y1 = control.step_response(flt, T=tt)
t1,y1 = control.step_response(700*plant, T=tt)
t3,y3 = control.step_response(700*plant*flt, T=tt)
c_p = control.feedback(ctrl*plant)
t2,y2 = control.step_response(200*c_p, T=tt)
plt.figure('Controller')
plt.plot(t1, y1, label='Plant')
plt.plot(t3, y3, label='Plant*Filter')
plt.plot(t2, y2, label='K*G')
plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
plt.ylim((0,400))
plt.legend()
plt.grid()

plt.figure('Bode')
control.bode_plot((plant, ctrl*plant), dB=True)
plt.legend(('plant','G*plant'))

plt.figure('Simulation')
plt.plot(tt, temp-1100, label='temp')
t3,y3,x3 = control.forced_response(c_p, T=tt, U=temp-1100)
plt.plot(t3, y3, label='U*K*G')
plt.plot(tt, temp_filtered_2-1100, label='temp_filtered_2')
t4,y4,x4 = control.forced_response(c_p, T=tt, U=temp_filtered_2-1100)
plt.plot(t4, y4, label='filt(U)*K*G')
plt.legend()
plt.grid()


ctrl2 = 0.2 * control.tf((1,43), (1,0)) * control.tf((1,7), (1/10,1))
c_p2 = control.feedback(ctrl2*plant)
t2,y2 = control.step_response(200*c_p2, T=tt)
plt.figure('Controller 2')
plt.plot(t1, y1, label='Plant')
plt.plot(t2, y2, label='K2*G')
plt.plot(tt, temp-1100, label='temp')
plt.legend()
plt.grid()

plt.figure('Bode 2')
control.bode_plot((plant, ctrl2, ctrl2*plant), dB=True)
plt.legend(('plant','G2','G2*plant'))

tsp = 200
sim_wfs = np.zeros(tt.size)
sim_temp = np.zeros(tt.size)
ctrl_state = np.zeros(2)
plant_state = 0
for i in range(1,tt.size):
    tsp = temp_filtered_2[i] - 1100
    t1,y1,x1 = control.forced_response(ctrl2, T=tt[i-1:i+1], U=tsp-sim_temp[i-1], X0=ctrl_state)
    sim_wfs[i] = y1[-1]
    ctrl_state = x1[:,1]
    t2, y2, x2 = control.forced_response(plant, T=tt[i-1:i+1], U=sim_wfs[i], X0=plant_state)
    sim_temp[i] = y2[-1]
    plant_state = x2[1]

plt.figure('Simulation 3')
plt.plot(tt, temp_filtered_2-1100, label='temp')
plt.plot(tt, sim_wfs, label='sim_wfs')
plt.plot(tt, sim_temp, label='sim_temp')
plt.legend()
plt.grid()


# discrete time simulation
plant_d = control.tf(.3, (1/10, 1), .001)
plant_d = control.matlab.c2d(plant, .001)
ctrl2_d = control.matlab.c2d(ctrl2, .001)

t,y=control.step_response(plant_d,T=(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
plt.plot(t,y)

import weld1.filter1
fplant = weld1.filter1.Filter1((0,0.002985),(1., -0.99004983))
fctrl = weld1.filter1.Filter1(ctrl2_d.num[0][0],ctrl2_d.den[0][0])
b, a = signal.butter(3, 2/500, 'low')
ftemp = weld1.filter1.Filter1(b, a)
ff_temp = np.zeros(temp.size)
for i in range(temp.size):
    ff_temp[i] = ftemp.step(temp[i])
plt.figure('Check filter')
plt.plot(tt, temp, label='temp')
plt.plot(tt, ff_temp, label='ff_temp')
plt.legend()
plt.grid()


tsp = 200
sim_wfs = np.zeros(tt.size)
sim_temp = np.zeros(tt.size)
fplant = weld1.filter1.Filter1((0,0.002985),(1., -0.99004983))
ctrl2 = 0.2 * control.tf((1,43), (1,0)) * control.tf((1,7), (1/10,1))
# ctrl2 = 2 * control.tf((1/1,1), (1,0)) * control.tf((0/1,1), (0/2,1))
ctrl2 = 0.25 * control.tf((1,20), (1,0)) * control.tf((1,10), (1/15,1))
ctrl2 = control.tf(3, 1) + control.tf(2, (1,0)) + control.tf((0,0), 1)
ctrl2 = control.tf(7.5, 1) + control.tf(50, (1,0)) + control.tf((0.25,0), (0,1))
ctrl2_d = mt.c2d(ctrl2, .001)
fctrl = weld1.filter1.Filter1(ctrl2_d.num[0][0],ctrl2_d.den[0][0])
b, a = signal.butter(3, 100/500, 'low')
ftemp = weld1.filter1.Filter1(b, a)
for i in range(1,tt.size):
    tsp = max(temp[i] - 1100, 0)
    # err = tsp - ftemp.step(sim_temp[i-1])
    err = tsp - sim_temp[i - 1]  # ftemp.step(sim_temp[i-1])
    sim_wfs[i] = fctrl.step(err)
    sim_temp[i] = fplant.step(sim_wfs[i])

plt.figure('Simulation 10')
plt.plot(tt, temp-1100, label='temp')
plt.plot(tt, sim_wfs, label='sim_wfs')
plt.plot(tt, sim_temp, label='sim_temp')
plt.legend()
plt.grid()

plt.ylim((-400,1200))


plant = control.tf(.3, (1/10, 1))
b, a = signal.butter(3, 2/500, 'low')
flt = control.tf(b, a, .001)
# t1,y1 = control.step_response(flt, T=tt)
t1,y1 = control.step_response(700*plant_d, T=mt.linspace(0,12,1201))
t3,y3 = control.step_response(700*plant_d*flt, T=mt.linspace(0,12,1201))
c_p = control.feedback(ctrl2_d*plant_d*flt)
t2,y2 = control.step_response(200*c_p, T=mt.linspace(0,12,1201))
plt.figure('Controller')
plt.plot(t1, y1, label='Plant')
plt.plot(t3, y3, label='Plant*Filter')
plt.plot(t2, y2, label='K*G')
# plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
plt.ylim((0,400))
plt.legend()
plt.grid()


Kp=3
Ki=2
Kd=0
ctrl_d = 1 * mt.c2d(control.tf(Kp, 1) + control.tf(Ki, (1,0)) + control.tf((Kd,0), 1), .001)
ctrl2 = 0.5 * control.tf((1,20), (1,0)) * control.tf((1,10), (1/15,1))
# ctrl2 = 2 * control.tf((1/1,1), (1,0)) * control.tf((0/1,1), (0/2,1))
ctrl_d = mt.c2d(ctrl2, .001)
plt.figure('Bode')
control.bode_plot((ctrl_d,plant_d,ctrl_d*plant_d), dB=True)
plt.legend(('ctrl','plant','ctrl_d*plant_d'))

plt.figure('Step')
t1,y1 = control.step_response(plant_d)
t3,y3 = control.step_response(plant_d*flt)
c_p = control.feedback(ctrl_d*plant_d*flt)
t2,y2 = control.step_response(c_p)
c_p2 = control.feedback(ctrl_d*plant_d)
t4,y4 = control.step_response(c_p2)
# plt.plot(t1, y1, label='Plant')
# plt.plot(t3, y3, label='Plant*Filter')
plt.plot(t2, y2, label='K*G*F')
plt.plot(t4, y4, label='K*G')
# plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
# plt.ylim((0,400))
plt.legend()
plt.grid()

Kp=7.5
Ki=50
Kd=0.25
ctrl = 1*(control.tf(Kp, 1) + control.tf(Ki, (1,0)) + control.tf((Kd,0), 1))
ctrl2 = 0.25 * control.tf((1,20), (1,0)) * control.tf((1,10), (1/15,1))
ctrl2 = 1*(control.tf(Kp, 1) + control.tf(Ki, (1,0)) + control.tf((Kd,0), (1/100,1)))
plt.figure('Bode comp')
control.bode_plot((ctrl, ctrl2), dB=True)
plt.legend(('ctrl','ctrl2'))

plt.figure('Comp ctrls')
t1,y1 = control.step_response(ctrl*plant)
t2,y2 = control.step_response(ctrl2*plant)
plt.plot(t2, y2, label='ft')
plt.plot(t1, y1, label='pid')
plt.legend()
plt.grid()
plt.figure('Bode comp')
control.bode_plot((ctrl*plant,ctrl2*plant,(control.tf(1, 1) + control.tf(1, (1,0)) + control.tf((.1,0), 1))*plant), dB=True)
plt.legend(('pid','ft','pid*plant','ft*plant'))


###################### PLC algrithm simulation

tsp = 200
sim_wfs = np.zeros(tt.size)
sim_temp = np.zeros(tt.size)
fplant = weld1.filter1.Filter1((0,0.002985),(1., -0.99004983))
ctrl2 = 0.2 * control.tf((1,43), (1,0)) * control.tf((1,7), (1/10,1))
# ctrl2 = 2 * control.tf((1/1,1), (1,0)) * control.tf((0/1,1), (0/2,1))
ctrl2 = 0.25 * control.tf((1,20), (1,0)) * control.tf((1,10), (1/15,1))
ctrl2 = control.tf(3, 1) + control.tf(2, (1,0)) + control.tf((0,0), 1)
ctrl2 = control.tf(7.5, 1) + control.tf(50, (1,0)) + control.tf((0.25,0), (0,1))
ctrl2_d = mt.c2d(ctrl2, .001)
fctrl = weld1.filter1.Filter1(ctrl2_d.num[0][0],ctrl2_d.den[0][0])


Kp=7.5
Ki=10
Kd=0.5
Fd=.1
int_error = 0
# last_error = 0
d_state = 0
plant = control.tf(.3, (1/10, 1))
plant = control.tf(.3, (1/5, 1))
plant_d = control.matlab.c2d(plant, .001)
fplant = weld1.filter1.Filter1((0,plant_d.num[0][0][0]),plant_d.den[0][0])
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

plt.figure('Plc Simulation 6')
plt.plot(tt, temp-1100, label='temp')
plt.plot(tt, sim_wfs, label='sim_wfs')
plt.plot(tt, sim_temp, label='sim_temp')
plt.legend()
plt.grid()
plt.ylim((0,1200))


plt.figure('Seguimiento consigna')
plt.title('Seguimiento consigna')
plt.plot(tt, temp, label='consigna')
# plt.plot(tt, sim_wfs, label='sim_wfs')
plt.plot(tt, sim_temp+1100, label='plant2')
plt.legend()
plt.grid()
plt.ylim((1000,1500))
plt.xlabel('Time (s)')
plt.ylabel('Temp (C)')