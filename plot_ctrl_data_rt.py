import matplotlib.pyplot as plt
import numpy as np
import paho.mqtt.client as mqtt
import time


ctrl_data = None

def on_message(mqttc, obj, msg):
    global ctrl_data
    # print(msg.topic + " " + str(msg.qos) + " " + str(msg.payload))
    if ctrl_data is None:
        ctrl_data = np.fromstring(msg.payload)
    else:
        ctrl_data = np.vstack((ctrl_data, np.fromstring(msg.payload)))
    # print('on_message() ctrl_data.shape', ctrl_data.shape)

mqtt_client = mqtt.Client()
mqtt_client.on_message = on_message
mqtt_client.connect('172.24.1.45')
mqtt_client.subscribe("ctrl3d", 0)
mqtt_client.loop_start()


fig = plt.figure('Control rt')
labels_unit_scale=(('Power','(w)',3000),('Speed','(m/s)',0.01),('T1','(ºC)',2000),('T2','(ºC)',1000))
axs,line,n=[],[],len(labels_unit_scale)
for i in range(n):
    ax = fig.add_subplot(100 * n + 11 + i)
    lin, = ax.plot([], [], label=labels_unit_scale[i][0])
    ax.set_ylabel(labels_unit_scale[i][1])
    ax.set_ylim((0, labels_unit_scale[i][2]))
    ax.legend()
    ax.grid()
    axs.append(ax)
    line.append(lin)
axs[0].set_title('Control temperature')
axs[3].set_xlabel('(s)')
fig.show(0)

while True:
    cnt = 1
    if ctrl_data is not None and ctrl_data.shape[0] > cnt:
        cnt = ctrl_data.shape[0]
        # print('while- ctrl_data.shape', ctrl_data.shape)
        for i in range(n):
            line[i].set_xdata(ctrl_data[:,0])
            line[i].set_ydata(ctrl_data[:, i+1])
            axs[i].set_xlim((ctrl_data[0,0],ctrl_data[-1,0]))
        axs[0].set_title('Control temperature t='+str(round(ctrl_data[-1,0],2)))
        fig.canvas.draw()
        fig.canvas.flush_events()
    else:
        time.sleep(.5)