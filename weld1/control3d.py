import matplotlib.pyplot as plt
import numpy as np
import control

import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics()

plant = fem2ss.get_ss((.025,.020,0))
t1,y1 = control.step_response(2500*plant)
plt.figure('Weld model')
plt.plot(t1, y1, label='Plant')
plt.xlim((0,12))