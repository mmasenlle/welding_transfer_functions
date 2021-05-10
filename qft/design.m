load('ss_vars_reverse');
sys_pert = ss(A,[B Bv T0r],[C1;C2],0);
rsysp = balred(sys_pert, 100);
systf = tf(rsysp);

M1 = balred(systf(:,3), 2);
impulse(systf(:,3), M1)

s=tf('s');
M=[100*(s^2 - 1.881*s + 1.417)/((s+0.4704)*(s+0.3255));
   400*(s+3.198)*(s+0.512)/((s+0.4704)*(s+0.3255))];
impulse(systf(:,3), M),legend('D(s)','M(s)')

P=systf(:,1:2);
D=systf(:,3);

U=P^-1 * (M-D); % NaN

U1=P(1,1)^1 * (M(1)-D(1));
U2=P(2,2)^1 * (M(2)-D(2));
U=[U1;U2];
impulse(U)