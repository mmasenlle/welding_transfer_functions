import matplotlib.pyplot as plt
import numpy as np

cd = np.load('ctrl_qft_manu3.npy')
# cd = np.load('data/ctrl_data_pert.npy')
fig = plt.figure('Control 1')
labels_unit_scale=(('Power','(w)',15000,0),('Speed','(m/s)',0.02,-0.02),('T1','(ºC)',2000,0),('T2','(ºC)',2000,0))
axs,line,n=[],[],len(labels_unit_scale)
for i in range(n):
    ax = fig.add_subplot(100 * n + 11 + i)
    lin, = ax.plot(cd[:,0], cd[:,i+1], label=labels_unit_scale[i][0])
    ax.set_ylabel(labels_unit_scale[i][1])
    ax.set_ylim((labels_unit_scale[i][3], labels_unit_scale[i][2]))
    ax.legend()
    ax.grid()
    axs.append(ax)
    line.append(lin)
axs[0].set_title('Control temperature')
axs[3].set_xlabel('(s)')

fig = plt.figure('Control 100')
plt.plot(cd[:,0], cd[:,1]/5.0, label='Power (w)')
plt.plot(cd[:,0], cd[:,2]*3e5, label='Speed (um/s/5)')
plt.plot(cd[:,0], cd[:,3], label='T1 (ºC)')
plt.plot(cd[:,0], cd[:,4], label='T2 (ºC)')
plt.title('Control change direction')
plt.xlabel('(s)')
plt.grid()
plt.legend()

fig = plt.figure('Actions')
plt.plot(cd[:,0], cd[:,1]/5.0, label='Power (w)')
plt.plot(cd[:,0], cd[:,2]*3e5, label='Speed (um/s/5)')
plt.plot(cd[:,0], cd[:,3], label='T1 (ºC)')
plt.plot(cd[:,0], cd[:,4], label='T2 (ºC)')
plt.title('Actions')
plt.xlabel('(s)')
plt.grid()
plt.legend()

fig = plt.figure('Filter actions')
plt.plot(cd[:,0], cd[:,5], label='Power (w)')
plt.plot(cd[:,0], cd[:,6]*3e5, label='Speed (um/s/5)')
plt.title('Control filter actions')
plt.xlabel('(s)')
plt.grid()
plt.legend()

cd = np.load('ctrl_qft_manu4.npy')
od = np.load('test_open/open_loop_data.npy')
fig = plt.figure('Open vs Control 1')
# plt.plot(cd[:,0], cd[:,1]/5.0, label='Power (w)')
# plt.plot(cd[:,0], cd[:,2]*3e5, label='Speed (um/s/5)')
plt.plot(cd[:,0], cd[:,3], label='Ctrl T1 (ºC)')
# plt.plot(cd[:,0], cd[:,4], label='Ctrl T2 (ºC)')
plt.plot(od[:,0], od[:,3], label='Open T1 (ºC)')
# plt.plot(od[:,0], od[:,4], label='Open T2 (ºC)')
plt.title('Open vs Control 1')
plt.xlabel('(s)')
plt.grid()
plt.legend()
fig = plt.figure('Open vs Control 2')
plt.plot(cd[:,0], cd[:,4], label='Ctrl T2 (ºC)')
plt.plot(od[:,0], od[:,4], label='Open T2 (ºC)')
plt.title('Open vs Control 2')
plt.xlabel('(s)')
plt.grid()
plt.legend()