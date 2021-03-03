
load('ss_vars_isotherms')
sys1200 = [ss(A,B,C_1200x,0) ss(A,B,C_1200y,0)
           ss(A,Bv,C_1200x,0) ss(A,Bv,C_1200y,0)];

plant1=sys1200(1,1);
[p1,t1]=step(plant1*7500);

plant2=sys1200(2,2);
[p2,t2]=step(plant2*-.005);

plot(t1,p1,t2,p2*.5,[0 8],[1200 1200], '--k')

resp1=tf([1200],[1/1 1],'InputDelay',2.0)
[r1,tr1]=step(resp1);
resp2=tf([1200],[1/3 1])
[r2,tr2]=step(resp2);

rtf = -0.0083424 * (s-3.051) * (s^2 - 4.368*s + 10.91) / ((s+0.9005) * (s+1) * (s^2 + 1.742*s + 1.842));
[r3,tr3]=step(rtf*7500); zpk(rtf*7500)

plot(t1,p1,tr1,r1,tr3,r3,t2,p2*.5,tr2,r2,[0 8],[1200 1200], '--k')
legend('step1','respuesta 1a','respuesta 1b','step2','respuesta 2')

hold on
plot(tr1,r1,'g')
plot(tr2,r2,'y')

action1=resp1*1/plant1;
action2=resp2*1/plant2;
figure
step(action1)
figure
step(action2)

%%%%
% Order reduction
rtf1 = tf(balred(plant1, 3));
rtf2 = tf(balred(plant1, 2));zpk(rtf2)
bode(plant1, rtf1)
step(rtf1*7500)
zpk(rtf1)
s=tf('s');
rtf = -0.0083424 * (s-3.051) * (s^2 - 4.368*s + 10.91) / ((s+0.9005) * (s^2 + 1.742*s + 1.842));
rtf = -0.0083424 * (s-3.051) * (s^2 - 4.368*s + 10.91) / ((s+0.9005) * (s+1) * (s^2 + 1.742*s + 1.842));
rtf = 0.036517 * (s^2 - 3.115*s + 3.476) / ((s+1) * (s^2 + 1.386*s + 0.7583));
step(rtf1, rtf)
step(plant1, rtf)
step(plant1, rtf2, rtf)

act1=plant1*7500/plant1;
act1=rtf*7500/plant1;
[y,t]=step(act1,10);
axis([1 10 0 100000])
