clear all; clc; close all;

load('ss_vars_2d_v001_q100000');
v=0.001;
q=100000;

% load('ss_vars_2d_v002_q100000');
% v=0.002;
% q=100000;

% load('ss_vars_2d_v003_q100000');
% v=0.003;
% q=100000;

% load('ss_vars_2d_v002_q50000');
% v=0.002;
% q=50000;

% load('ss_vars_2d_v002_q200000');
% v=0.002;
% q=200000;

w = logspace(-3,2,100);

k=24.0;
rho=7925;
cp=460;
a=k/(rho*cp);
h=1000;
b=2*h/(rho*cp);
s = 1j*w;


R = @(x,y)  sqrt(x^2+y^2);
D = @(s) sqrt(4*a*s + 4*a*b + v^2);

chi = @(x,y,s) v*besselk(0,R(x,y)*D(s)/(2*a)) - (-x/R(x,y))*(2*abs(v)-(v^2)./D(s)).*besselk(1,R(x,y)*D(s)/(2*a));
ft_vel = @(x,y) q*exp(-v*x/(2*a))./(4*pi*k*a*s).*(chi(x,y,s) - chi(x,y,0));


% comparativa ft ss +/-
figure
bode(frd(ft_vel(0,0.005), w), ss(A,Bv,C_00__05,0), frd(ft_vel(0,-0.005), w), ss(A,Bv,C_00_05,0), w)
title('Bode T/Vel y=5,-5')
legend('ft+','ss+','ft-','ss-')

% puntos longitudinales
figure
bode(frd(ft_vel(-0.005,0), w), frd(ft_vel(-0.01,0), w), frd(ft_vel(-0.02,0), w), ss(A,Bv,C_05_00,0),ss(A,Bv,C_10_00,0),ss(A,Bv,C_20_00,0),w)
title('Bode T/Vel x=-5,-10,-20')
legend('ft -5,0','ft -10,0','ft -20,0','ss -5,0','ss -10,0','ss -20,0')
figure
bode(frd(ft_vel(0.003,0), w), frd(ft_vel(0.005,0), w), frd(ft_vel(0.007,0), w), ss(A,Bv,C__03_00,0),ss(A,Bv,C__05_00,0),ss(A,Bv,C__07_00,0),w)
title('Bode T/Vel x=3,5,7')
legend('ft 3,0','ft 5,0','ft 7,0','ss 3,0','ss 5,0','ss 7,0')

% puntos laterales
figure
bode(frd(ft_vel(0,0.005), w), frd(ft_vel(0,0.007), w), frd(ft_vel(0,0.01), w), ss(A,Bv,C_00_05,0),ss(A,Bv,C_00_07,0),ss(A,Bv,C_00_10,0),w)
title('Bode T/Vel  y=5,7,10')
legend('0,5','0,7','0,10')
