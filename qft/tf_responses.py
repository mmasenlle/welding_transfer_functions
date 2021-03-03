import matplotlib.pyplot as plt
import numpy as np
import control


plant = control.tf((1200,),(1/5,1))
t1,y1 = control.step_response(plant)
plt.figure('Responses')
plt.plot(t1, y1, label='Response1')
# plt.xlim((0,12))
