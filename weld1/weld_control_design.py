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
plt.plot(tt,temp, label='Temp (C)')
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
plt.figure('Plant')
plt.plot(t, 700*y, label='plant')
plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
plt.ylim((0,400))
plt.legend()


import control
Kp=10
Ki=60
Kd=0.2
ctrl = control.tf(Kp, 1) + control.tf(Ki, (1,0)) + control.tf((Kd,0), 1)
plant = control.tf(.3, (1/10, 1))
t1,y1 = control.step_response(700*plant, T=tt)
c_p = control.feedback(ctrl*plant)
t2,y2 = control.step_response(200*c_p, T=tt)
plt.figure('Controller')
plt.plot(t1, y1, label='Plant')
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

