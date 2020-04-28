import numpy as np
import matplotlib.pyplot as plt
import control

k=24.0
rho=7925
cp=460
v=0.002
a=k/(rho*cp)
x=.01
y=0
z = 0
q=500
dl=0.001

w = np.logspace(-3,1,100)
s = 1j*w

ft_pot = lambda x,y,z: 1/(rho*cp)*np.exp(-(v*x+np.sqrt((x**2+y**2+z**2)*(4*a*s+v**2)))/(2*a))/(2*np.pi*a*np.sqrt(x**2+y**2+z**2))
plt.figure('Bode T_i/power')
control.bode_plot([control.frd(ft_pot(dl,.0,0), w),control.frd(ft_pot(-dl,.0,0), w),control.frd(ft_pot(.0,dl,0), w),control.frd(ft_pot(.0,0,dl), w),\
                   control.frd(ft_pot(.0,-dl,0), w),control.frd(ft_pot(.0,0,-dl), w)], w, dB=True)
plt.legend(['x='+str(dl)+'','x=-'+str(dl)+'','y='+str(dl)+'','z='+str(dl)+'','y=-'+str(dl)+'','z=-'+str(dl)+''])
plt.figure('Bode T_i/speed')
ft_vel = lambda x,y,z: -q*(v*np.sqrt(x**2 + y**2 + z**2) + x*np.sqrt(4*a*s + v**2))*np.exp(-(v*x + np.sqrt(4*a*s + v**2)*np.sqrt(x**2 + y**2 + z**2))/(2*a))/(4*np.pi*a**2*np.sqrt(4*a*s + v**2)*np.sqrt(x**2 + y**2 + z**2))
control.bode_plot([control.frd(ft_vel(dl,.0,0), w),control.frd(ft_vel(-dl,.0,0), w),control.frd(ft_vel(.0,dl,0), w),control.frd(ft_vel(.0,0,dl), w),\
                   control.frd(ft_vel(.0,-dl,0), w),control.frd(ft_vel(.0,0,-dl), w)], w, dB=True)
plt.legend(['x='+str(dl)+'','x=-'+str(dl)+'','y='+str(dl)+'','z='+str(dl)+'','y=-'+str(dl)+'','z=-'+str(dl)+''])ç


plt.figure('Bodes')
control.bode_plot([control.frd(ft_pot(.0,dl,0), w),control.frd(ft_vel(.0,dl,0), w),\
                   control.frd(ft_pot(.0,-dl,0), w),control.frd(ft_vel(.0,-dl,0), w)], w, dB=True)
plt.legend(['pot y='+str(dl)+'','vel y='+str(dl)+'','pot y=-'+str(dl)+'','vel y=-'+str(dl)])