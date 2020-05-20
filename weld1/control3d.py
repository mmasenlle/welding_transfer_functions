import matplotlib.pyplot as plt
import numpy as np
import control

import fem3d_2ss
fem2ss = fem3d_2ss.Fem3d_fenics()

plant = fem2ss.get_ss((.02,.025,0))
t1,y1 = control.step_response(1500*plant)
plt.figure('Weld model')
plt.plot(t1, y1, label='Plant')
plt.xlim((0,12))

plt.figure('Weld model plants')
# pts = ((.02,.02,0),(.021,.02,0),(.022,.02,0),(.023,.02,0),(.024,.02,0),(.025,.02,0),(.02,.021,0),(.02,.022,0),
#           (.02,.023,0),(.024,.02,0),(.025,.02,0),(.021,.021,0),(.022,.022,0),(.023,.023,0)):
# pts = ((.02, .02, 0), (.02, .021, 0), (.02, .022, 0), (.02, .023, 0), (.02, .024, 0), (.02, .025, 0), (.02, .026, 0)):
# pts = ((.02, .02, 0), (.021, .021, 0), (.022, .022, 0), (.023, .023, 0), (.024, .024, 0), (.025, .025, 0), (.026, .026, 0)):
pts = ((.02, .02, 0), (.021, .02, 0), (.02, .021, 0), (.02, .02, 0.001), (.019, .02, 0), (.02, .019, 0.0), (.021, .021, 0.001))
pts = ((.025, .02, 0), (.02, .025, 0), (.03, .02, 0.0), (.02, .03, 0), (.04, .02, 0.0), (.025, .025, 0.0))
for p in pts:
    plant = fem2ss.get_ss(p)
    t1, y1 = control.step_response(1500 * plant)
    plt.plot(t1, y1, label='Plant '+str(p))
plt.xlim((0,12))
plt.legend()



ctrl2 = control.tf2ss(1 * control.tf((1/1,1), (1,0)) * control.tf((0/1,1), (0/2,1)))
# t1,y1 = control.step_response(flt, T=tt)
t1,y1 = control.step_response(plant)
c_p = control.feedback(ctrl2*plant)
t2,y2 = control.step_response(c_p)
plt.figure('Bode')
control.bode_plot((ctrl2,plant,ctrl2*plant), dB=True)
plt.legend(('ctrl','plant','ctrl*plant'))
plt.figure('Controller 2')
plt.plot(t1, y1, label='Plant')
plt.plot(t2, y2, label='K*G')
# plt.plot(tt, temp-1100, label='temp')
# plt.plot(tt, temp_filtered_2, label='temp_filtered_2')
# plt.ylim((0,400))
plt.legend()
plt.grid()