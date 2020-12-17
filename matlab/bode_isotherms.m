clear all; clc; close all;

bode_opt = bodeoptions;
set(bode_opt,'PhaseMatching','on')

load('ss_vars_isotherm')


figure
bode(ss(A,B,C_1200x,0),ss(A,B,C_1000x,0),ss(A,B,C_800x,0),ss(A,B,C_600x,0))
title('Bode T1/Pot isotherms')
legend('1200','1000','800','600')

figure
bode(ss(A,Bv,C_1200y,0),ss(A,Bv,C_1000y,0),ss(A,Bv,C_800y,0),ss(A,Bv,C_600y,0))
title('Bode T2/Vel isotherms')
legend('1200','1000','800','600')