load('ss_vars_reverse');
sys = ss(A,[B Bv],[C1;C2],0);

figure
subplot(4,1,1),plot(od(:,1),od(:,2));title('Non lineal open loop Power')
subplot(4,1,2),plot(od(:,1),od(:,3));title('Non lineal open loop Speed')
subplot(4,1,3),plot(od(:,1),od(:,4));title('Non lineal open loop T1')
subplot(4,1,4),plot(od(:,1),od(:,5));title('Non lineal open loop T2')

t=od(:,1);
u=[od(:,2) od(:,3)*0];
y = lsim(sys,u,t,T0r);
figure
subplot(4,1,1),plot(t,u(:,1));title('Lineal open loop Power')
subplot(4,1,2),plot(t,u(:,2));title('Lineal open loop Speed')
subplot(4,1,3),plot(t,y(:,1));title('Lineal open loop T1')
subplot(4,1,4),plot(t,y(:,2));title('Lineal open loop T2')

% comparation
figure
subplot(2,1,1),plot(t,od(:,4),t,y(:,1));title('T1'); legend('FEM','ss')
subplot(2,1,2),plot(t,od(:,5),t,y(:,2));title('T2'); legend('FEM','ss')

% s=tf('s')
% sia = s*eye(length(sys.A)) - sys.A;
% isia = inv(s*eye(length(sys.A)) - sys.A);
% D=sys.C * inv(s*eye(length(sys.A)) - sys.A) * T0r;


load('ss_vars_reverse');
sys_pert = ss(A,[B Bv T0r],[C1;C2],0);
t=od(:,1);
u3=zeros(length(t),1); u3(1)=1/(t(2)-t(1));
u=[od(:,2) od(:,3)*0 u3];
y = lsim(sys_pert,u,t);
figure
subplot(2,1,1),plot(t,od(:,4),t,y(:,1));title('T1'); legend('FEM','ss')
subplot(2,1,2),plot(t,od(:,5),t,y(:,2));title('T2'); legend('FEM','ss')

% reduccion de sistema
rsysp = balred(sys_pert, 100);
y = lsim(rsysp,u,t);
figure
subplot(2,1,1),plot(t,od(:,4),t,y(:,1));title('T1'); legend('FEM','ss red')
subplot(2,1,2),plot(t,od(:,5),t,y(:,2));title('T2'); legend('FEM','ss red')

% modo ft
systf = tf(rsysp);
y = lsim(systf,u,t);
figure
subplot(2,1,1),plot(t,od(:,4),t,y(:,1));title('T1'); legend('FEM','tf')
subplot(2,1,2),plot(t,od(:,5),t,y(:,2));title('T2'); legend('FEM','tf')

figure; bode(systf)
impulse(systf(:,3))
