
% load('ss_vars_isotherms')
% sys1200 = [ss(A,B,C_1200x,0) ss(A,Bv,C_1200x,0)
%            ss(A,B,C_1200y,0) ss(A,Bv,C_1200y,0)];

% rsys = tf(balred(sys1200, 3));
% save('reduced_plant.mat', 'rsys')

load('reduced_plant')
% step(sys1200, rsys);

s=tf('s');

rsysd = [rsys(1,1) 0; 0 rsys(2,2)];
step([7500 0; 0 -.0025]*rsysd)

Gc=[(s+1)/s 0; 0 (-5e-06*s-1.5e-05)/(3*s)];
%step([1200 0; 0 1200]*feedback(Gc*rsys, eye(2)))
%impulse([1200/s 0; 0 1200/s]*feedback(Gc*rsys, eye(2)))

resp=[1200/s 0; 0 1200/s]*feedback(Gc*rsysd, eye(2));
impulse(resp)

resp2=resp*[1;1]; %[resp(1,1); resp(2,2)];
impulse(resp2, 50)
title('Respuestas del sistema a referencias escalón')

figure
act=rsysd^-1 * resp2;
impulse(act, 50)
title('Acciones de control para referencias escalón :(')
% no sale muy bien ????

resp_defined = [
 7200 * (-0.0083424 * (s-3.051) * (s^2 - 4.368*s + 10.91)) / (s * (3*s+0.9005) * (3*s+1) * (s^2 + 1.742*s + 1.842));  
 1200/(s*(s/2 + 1))];
figure
subplot(2,1,1),impulse(7500/s * rsys(1,1), resp2(1), resp_defined(1))
title('Respuestas del sistema')
subplot(2,1,2),impulse(-.0025/s * rsys(2,2), resp2(2), resp_defined(2))
legend('Acción escalón', 'Seguimiento escalón', 'Definida')
title('')

figure
act2=rsys^-1 * resp_defined;
impulse(act2, 1)
title('Acciones de control para las respuestas deseadas :(')


% pruebas SISO
% impulse(7500/s * rsys(1,1) / rsys(1,1))
% impulse(resp(1,1) / rsys(1,1))
% impulse(resp_defined(1) / rsys(1,1))

% impulse(-.0025/s * rsys(2,2) / rsys(2,2))
% impulse(resp(2,2) / rsys(2,2))
% impulse(resp_defined(2) / rsys(1,1))

