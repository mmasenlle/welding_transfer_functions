
clear all; close all; clc;

% load('ss_plants2')
% for k = 1:size(AAA,3)
%     delta_T0 = TT0r(:,:,k)-TT0(:,:,k);
%     sys_red(:,:,k) = balred(ss(AAA(:,:,k),[BB(:,k) BBv(:,k) delta_T0],[C1;C2],0),40);
% end
% sys_red.u = {'P','V','delta'};
% sys_red.y = {'T1','T2'};

load('plantas')

figure(1); clf; step(sys_red(:,3,:))

P = sys_red(:,1:2,:);
D = sys_red(:,3,:);
P_som = inv(P);
nplants = size(P,3);

% Especificaciones 
s = tf('s');
% b11_spec = -3600*(s/2.5+1)/(s/1+1)/(s/1+1);
% b21_spec = 4000*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);

% relajacion de especificaciones ?
b11_spec = -3000*(s/2.5+1)/(s/1+1)/(s/1+1);
b21_spec = 4200*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);

% seleccion de frecuencias de diseño
figure(1); clf; bode(sys_red(1,1,:),b11_spec)
figure(1); clf; bode(sys_red(2,2,:),b21_spec)

w = [0.01 0.1 1 3 5 10 30 100];
for ki = 1:nplants
    kf = 0;
    for fr = w
        kf = kf+1;
        P_som_aux = freqresp(P_som(:,:,ki),fr);
        p11_som(ki,kf) = P_som_aux(1,1);
        p12_som(ki,kf) = P_som_aux(1,2);
        p21_som(ki,kf) = P_som_aux(2,1);
        p22_som(ki,kf) = P_som_aux(2,2);
        D_aux = freqresp(D(:,:,ki),fr);
        d11(ki,kf) = D_aux(1,1);
        d21(ki,kf) = D_aux(2,1);
    end
end

Bs = [b11_spec; b21_spec];
% figure;bodemag(Bs,D) %Para ver si las especificaciones son alcanzables

w1=logspace(-3,2,2000);
phs=0:-0.1:-360;

nom = 14; % planta nominal benchmark
P110 = 1/P_som(1,1,nom);
% plottmpl(w,1/P_som(1,1),nom);

b11_new = abs(squeeze(freqresp(Bs(1,1),w))');
b21_new = abs(squeeze(freqresp(Bs(2,1),w))');

% Factores
mult11_new = [1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00];
mult21_new = [1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00];

b11 = mult11_new.*b11_new;
b21 = mult21_new.*b21_new;

% x11 = [0.25 0.25 0.5 0.5 0.5 0.5 0.5 0.5];
% x11 = [0.3  0.3  0.55 0.5 0.5  0.1  0.1 0.2];
% x11 = [0.47 0.48 0.21 0.22 0.23 0.11 0.07];
x11 = [0.472  0.5  0.55 0.3 0.201  0.1  0.1 0.2];
x21 = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];

bc11 = b11.*x11;
bc21 = b21.*x21; %dudas

b11_aux = repmat(b11,nplants,1); %se hacen copias (en filas) para poder operar: una para cada planta dentro de la incertidumbre
b21_aux = repmat(b21,nplants,1); 

%% DiseÃ±o del lazo 1
PM = 30; %Margen de fase robusto objetivo
Wstab1 = 0.5/(cos(pi*(180-PM)/(2*180)));
B1_stab = sisobnds(1,w,Wstab1,1/P_som(1,1),[],1,[],1,phs);

%las columnas son las frecuencias de diseÃ±o, las filas son la incertidumbre
aux11 = freqresp(P_som(1,1),w);
aux12 = freqresp(P_som(1,2),w);
p11_aux = reshape(aux11,size(aux11,3),size(aux11,4)).';
p12_aux = reshape(aux12,size(aux12,3),size(aux12,4)).';
d11e_p11 = abs(p12_aux).*b21_aux./p11_aux;
B11_coup = genbnds(10,w,bc11,d11e_p11,0,1,1/P_som(1,1),P110,phs);

% QD11_tf = minreal(-P_som(1,1)*D_mat(1,1)-P_som(1,2)*D_mat(2,1));
% B11_distff_tf = sisobnd14c(w,tf(0,1),b11-bc11,1/P_som(1,1),QD11_tf,nom,phs);
% B11_distff_ad_tf = adaptation(B11_distff_tf,14);

QD11_fr = -p11_som.*d11-p12_som.*d21;
B11_distff_fr = sisobnd14c(w,tf(0,1),b11-bc11,1./p11_som,QD11_fr,nom,phs);
B11_distff_ad_fr = adaptation(B11_distff_fr,14);

% plotbnds(B11_distff_ad_tf,[],phs);
% plotbnds(B11_distff_ad_fr,[],phs);

B11_distff_ad = B11_distff_ad_fr;
Bc11 = grpbnds(B1_stab,B11_coup);%,B11_distff_ad);
% Bs11 = sectbnds(Bc11);

% plotbnds(Bc11,[],phs);
% plotbnds(Bs11,[],phs);
% set(gcf,'units','points','position',[700 300 375 300])
% set(gca,'position',[0.1 0.1 0.85 0.85],'ylim',[-40 30]); grid off; legend off;

close all;
G110 = tf(1.8,[1 0]);
lpshape(w1,Bc11,P110,G110,phs);
% lpshape(w1,Bc11,1/p110_frd,G110,phs);


% % DiseÃ±o de f11 (FT)
% phs_F = 0:-0.01:-360;
% A11 = minreal(QD11_tf/P_som(1,1));
% B11 = minreal(-1/P_som(1,1));
% C11 = minreal(1+G110/P_som(1,1));
% D11 = 0;
% Bf11 = genbnds(10,w,b11-bc11,A11,B11,C11,D11,tf(1,1),phs_F);
% F110 = tf(-0.59*[1/0.03 1],[1/0.012 1])*tf(1,[1/0.3 1]); %Para Jor1
% lpshape(w1,Bf11,tf(1,1),F110,phs_F);

% DiseÃ±o de f11 (FR)
g11 = squeeze(freqresp(G110,w)).';
g11_fr = repmat(g11,nplants,1);
phs_F = 0:-0.01:-360;
A11 = QD11_fr./p11_som;
B11 = -1./p11_som;
C11 = 1+g11_fr./p11_som;
D11 = 0;
Bf11 = genbnds(10,w,b11-bc11,A11,B11,C11,D11,tf(1,1),phs_F);
% F110 = -97426.5/(s/1+1)/(s/0.23854+1); %
% F110 = -minreal(P_som(1,1,nom)*D(1,1,nom)+P_som(1,2,nom)*D(2,1,nom));
F110 = -7374.9064*(s/2.7093+1)*(s/-0.051239+1)*(s/0.30824+1)*(s^2/2.6225^2+2*s*0.40054/2.6225+1)*(s^2/13.966^2+2*s*0.71187/13.966+1)/(s/19.6814+1)/(s^2/0.20329^2+2*s*0.69779/0.20329+1)/(s^2/2.144^2+2*s*0.11508/2.144+1)/(s^2/7.5866^2+2*s*0.61737/7.5866+1);
lpshape(w1,Bf11,tf(1,1),F110,phs_F);

%% DiseÃ±o del lazo 2

f11 = squeeze(freqresp(F110,w)).';
f11_fr = repmat(f11,nplants,1);

P_som2(2,1,:) = minreal(P_som(2,1) - (P_som(2,1)*P_som(1,1))/(P_som(1,1)+G110));
P_som2(2,2,:) = minreal(P_som(2,2) - (P_som(2,1)*P_som(1,2))/(P_som(1,1)+G110));

p21_som2 = p21_som - p21_som.*p11_som./(p11_som+g11_fr);
p22_som2 = p22_som - p21_som.*p12_som./(p11_som+g11_fr);

f11_tf2 =  minreal(-(P_som(2,1)*F110)/(P_som(1,1)+G110));
f11_fr2 =  -p21_som.*f11_fr./(p11_som+g11_fr);

Wstab2 = 0.5/(cos(pi*(180-30)/(2*180)));
B2_stab = sisobnds(1,w,Wstab2,1/P_som2(2,2),[],1,[],1,phs);

QD21_previo = minreal(-P_som(2,1)*D_mat(1,1)-P_som(2,2)*D_mat(2,1));
QD21_tf = minreal(-P_som2(2,1)*D_mat(1,1)-P_som2(2,2)*D_mat(2,1)-f11_tf2);
bode(QD21_previo,QD21_tf)
% B21_distff_tf = sisobnd14c(w,tf(0,1),b21-bc21,1/P_som2(2,2),QD21_tf,nom,phs);
% B21_distff_ad_tf = adaptation(B21_distff_tf,14);

QD21_fr = -p21_som2.*d11-p22_som2.*d21-f11_fr2;
B21_distff_fr = sisobnd14c(w,tf(0,1),b21-bc21,1./p22_som2,QD21_fr,nom,phs);
B21_distff_ad_fr = adaptation(B21_distff_fr,14);

% plotbnds(B21_distff_ad_tf,[],phs);
% plotbnds(B21_distff_ad_fr,[],phs);

B21_distff_ad = B21_distff_ad_fr;
Bc22 = grpbnds(B2_stab,B21_distff_ad);
% Bs22 = sectbnds(Bc22);

% plotbnds(Bc22,[],phs);
% plotbnds(Bs22,[],phs);
% set(gcf,'units','points','position',[700 300 375 300])
% set(gca,'position',[0.1 0.1 0.85 0.85],'ylim',[-40 30]); grid off; legend off;

P220 = 1/P_som2(2,2,nom);

close all;
G220 = tf(-9.3e-6,[1 0]);
lpshape(w1,Bc22,P220,G220,phs);

% % DiseÃ±o de f21 (TF)
% phs_F=0:-0.01:-360;
% A21 = minreal(QD21_tf/P_som2(2,2));
% B21 = minreal(-1/P_som2(2,2));
% C21 = minreal(1+G220/P_som2(2,2));
% D21 = 0;
% Bf21 = genbnds(10,w,b21-bc21,A21,B21,C21,D21,tf(1,1),phs_F);
% F210 = tf(0.3,[1/0.017 1])*tf(1,[1/0.044 1]);
% lpshape(w1,Bf21,tf(1,1),F210,phs_F);

% DiseÃ±o de f11 (FR)
g22 = squeeze(freqresp(G220,w)).';
g22_fr = repmat(g22,nplants,1);
phs_F = 0:-0.01:-360;
A21 = QD21_fr./p22_som2;
B21 = -1./p22_som2;
C21 = 1+g22_fr./p22_som2;
D21 = 0;
Bf21 = genbnds(10,w,b21-bc21,A21,B21,C21,D21,tf(1,1),phs_F);
F210 = 0.0060744*(s/8.4199+1)/(s/1+1);
lpshape(w1,Bf21,tf(1,1),F210,phs_F);

%%
%%%%% ComprobaciÃ³n de resultados en la frecuencia %%%%%
close all;

Bs = [b11_spec; b21_spec];

G = [G110 0;0 G220];
F = [F110;F210];

ncomp = 100;
wcomp = logspace(-2,2,ncomp);

P_resp = freqresp(P(1:2,1:2),wcomp);
D_resp = freqresp(D,wcomp);
for kf = 1:ncomp
    G_resp = freqresp(G,wcomp(kf));
    F_resp = freqresp(F,wcomp(kf));
    for kp = 1:nplants
        L_resp(:,:,kf,kp) = P_resp(:,:,kf,kp)*G_resp;
        Tew_resp(:,:,kf,kp) = inv(eye(2)+L_resp(:,:,kf,kp))*(-D_resp(:,:,kf,kp)-P_resp(:,:,kf,kp)*F_resp);
    end
end

% E_frd = frd(E_resp,wcomp);
% figure; bodemag(E_frd,Bs,wcomp);

Tew_frd = frd(Tew_resp,wcomp);
figure; bodemag(Tew_frd,Bs,wcomp);

% hold on
% Tew = minreal(inv(eye(2)+P(1:2,1:2,:)*G)*(-D_mat-P(1:2,1:2,:)*F))
% bodemag(Tew)

% figure
% Tyw = minreal(inv(eye(2)+P(1:2,1:2,:)*G)*(D_mat+P(1:2,1:2,:)*F))
% step(Tyw)
% hold on


%% %%% ComprobaciÃ³n de resultados en el tiempo %%%%%

% G110 = tf(0.06,[1 0])*tf([1/0.025 1],[1/0.25 1]); %Jor1
% G220 = tf(0.06,[1 0])*tf([1/0.02 1],[1/1.2 1]); %Jorg1
% F110 = tf(-0.59*[1/0.03 1],[1/0.012 1])*tf(1,[1/0.3 1]); %Para Jor1
% F210 = tf(0.3*[1/0.075^2 2*0.25/0.075 1],[1/0.15^2 2*0.78/0.15 1])*tf(1,[1/0.3 1])*tf(1,[1/0.012 1]); %Para Jor1

[numG11,denG11]=tfdata(G110,'v');
[numG22,denG22]=tfdata(G220,'v');
[numF11,denF11]=tfdata(F110,'v');
[numF21,denF21]=tfdata(F210,'v');

% ParÃ¡metros de Johansson MP
A1 = 28; A2 = 32; A3 = 28; A4 = 32;
a1 = 0.071; a2 = 0.057; a3 = 0.071; a4 = 0.057; 
k1 = 3.33; k2 = 3.35; 
kc = 0.50; 
g = 981;
gamma1 = 0.7; gamma2 = 0.6;

amp = 0.5;

tsim = 300;
tref = linspace(0,tsim,1000);
[yb,tb] = step(amp*Bs,tref);

% close all
figure
set(gcf,'position',[400 300 300 400])

for caso = 1:nplants
    [numP11,denP11] = tfdata(P(1,1,caso),'v'); % Preparamos la simulaciÃ³n lineal
    [numP12,denP12] = tfdata(P(1,2,caso),'v');
    [numP21,denP21] = tfdata(P(2,1,caso),'v');
    [numP22,denP22] = tfdata(P(2,2,caso),'v');
    [numD11,denD11] = tfdata(P(1,3,caso),'v');
    [numD21,denD21] = tfdata(P(2,3,caso),'v');
    h10 = H_ss(1,caso); % Preparamos la simulaciÃ³n no lineal
    h20 = H_ss(2,caso);
    h30 = H_ss(3,caso);
    h40 = H_ss(4,caso);
    v10 = V_ss(1,caso);
    v20 = V_ss(2,caso);
    tsp1 = 0; amp1 = 0; tsp2 = 0; amp2 = 0; tsp3 = 0; amp3 = amp;
    sim('simulador_no_lineal')
    subplot(2,1,1);plot(t,y1-h10*kc,'color',[0.75 0.75 0.75]); hold on;
    subplot(2,1,2);plot(t,y2-h20*kc,'color',[0.75 0.75 0.75]); hold on;
%     sim simulador_lineal;
%     subplot(2,1,1);plot(t,y1,'g');
%     subplot(2,1,2);plot(t,y2,'g'); 
end

subplot(2,1,1); 
set(gca,'xlim',[0 tsim],'FontSize',9); %,'ylim',[-1 1]
plot(tb,yb(:,1,1),'--','color',[0 0 0],'linewidth',0.7);
plot(tb,-yb(:,1,1),'--','color',[0 0 0],'linewidth',0.7);

subplot(2,1,2); 
set(gca,'xlim',[0 tsim],'FontSize',9);%,'ylim',[-1 1]
plot(tb,yb(:,2,1),'--','color',[0 0 0],'linewidth',0.7);
plot(tb,-yb(:,2,1),'--','color',[0 0 0],'linewidth',0.7);

subplot(2,1,1); set(gca,'position',[0.11 0.55 0.85 0.43])
set(gca,'ytick',[-0.05 0 0.05])
subplot(2,1,2); set(gca,'position',[0.11 0.05 0.85 0.43])
set(gca,'ylim',[-0.02 0.02],'ytick',[-0.02 0 0.02])

%%
G110_jor = tf(0.06,[1 0])*tf([1/0.025 1],[1/0.25 1]); %Jor1
G220_jor = tf(0.06,[1 0])*tf([1/0.02 1],[1/1.2 1]); %Jorg1
G110_joh = tf(0.1,[1 0])*tf([1/0.0333 1],1); %Johansson
G220_joh = tf(0.0675,[1 0])*tf([1/0.0250 1],1); %Johansson
close all
bodemag(G110_jor,G110_joh)
close all
bodemag(G220_jor,G220_joh)













