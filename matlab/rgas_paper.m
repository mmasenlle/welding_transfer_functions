clear all; clc; close all;

%%%%%%%% Mallado irregular (1231 nodos)
load('ss_vars_irregular')


sys1 = [ss(A,B,C_10_00,0) ss(A,B,C_00_05,0)
        ss(A,Bv,C_10_00,0) ss(A,Bv,C_00_05,0)];


%% Design MF=30
w = 1; % Como no hay incertidumbre, basta una frecuencia (cualquiera)
W1 = 0.5/cos(150*pi/180/2); %MF=30
phs = -360:1:0;
w1=logspace(-3,2,1000);
s = tf('s');
% lazo 1
% B_stab = sisobnds(1,w,W1,ss(A,B,C_10_00,0),0,1,[],1,phs);
Gc_pot_10_00 = 2/s^1*(s/1.3+1);
% lpshape(w1,B_stab,ss(A,B,C_10_00,0),Gc_pot_10_00,phs)
% lazo 2
% B_stab=sisobnds(1,w,W1,ss(A,Bv,C_00_05,0),0,1,[],1,phs)
Gc_vel_00_05 = -0.00001/s^1*(s/4.453+1)/(s/10.2095+1);
% lpshape(w1,B_stab,ss(A,Bv,C_00_05,0),Gc_vel_00_05,phs)

%% Design MF=40
Gc_pot_10_00 = 1/s^1*(s/3+1);
Gc_vel_00_05 = -5e-06/s/(s/1+1);

figure; margin(Gc_pot_10_00*sys1(1,1));
figure; margin(Gc_vel_00_05*sys1(2,2));

%% Responses
sys1_diag = [sys1(1,1) 0
            0 sys1(2,2)];
Gc = [Gc_pot_10_00 0; 0 Gc_vel_00_05];
closed_loop_mimo = feedback(Gc * sys1, eye(2));
closed_loop_siso = feedback(Gc * sys1_diag, eye(2));
figure;
bodemag(closed_loop_mimo,closed_loop_siso);
title('Bode Plots of MIMO & SISO closed loops');
legend 'MIMO' 'SISO'

figure;
step(closed_loop_mimo,closed_loop_siso)
title('Steps of MIMO & SISO closed loops');
legend 'MIMO' 'SISO'

%% RGA analysis
% linsys=sys1d_4_1 
close all
w = logspace(-3,2,100);
Pf1 = freqresp(sys1,w);


for i=1:length(w)
    RGAw1(:,:,i)=Pf1(:,:,i).*inv(Pf1(:,:,i)).';
    RGAnum1(:,:,i) = sum(sum(abs(RGAw1(:,:,i) - eye(2))));
end

%squeeze(RGAnum)
AB=[0.5 0.5];
%figure(4); clf
figure
subplot(2,2,1); semilogx(w,abs(squeeze(RGAw1(1,1,:))))
     xlim([1e-3,1e1]);ylabel('RGA');%;xlabel('Frecuency (rad/s)')
    title('$$ RGAw(P,T_{1}) $$','interpreter','latex');grid;
%     legend('0.015','0.02','0.03','0.04','0.05');
    hold on; plot(AB,get(gca,'ylim'),'--y')  
    subplot(2,2,2); semilogx(w,abs(squeeze(RGAw1(1,2,:))))
     xlim([1e-3,1e1]);%xlabel('Frecuency (rad/s)');ylabel('RGA')
    title('$$ RGAw(P,T_{2}) $$','interpreter','latex');grid;
%     legend('x_2 = -0.02','x_2 = -0.03','x_2 = -0.04','x_2 = -0.05');
    hold on; plot(AB,get(gca,'ylim'),'--y')
    subplot(2,2,3); semilogx(w,abs(squeeze(RGAw1(2,1,:))))
     xlim([1e-3,1e1]);xlabel('Frecuency (rad/s)');ylabel('RGA')
    title('$$ RGAw(V,T_{1}) $$','interpreter','latex');grid;
%     legend('0.015','0.02','0.03','0.04','0.05');
    hold on; plot(AB,get(gca,'ylim'),'--y')
    subplot(2,2,4); semilogx(w,abs(squeeze(RGAw1(2,2,:))))
     xlim([1e-3,1e1]);xlabel('Frecuency (rad/s)');%ylabel('RGA')
    title('$$ RGAw(V,T_{2}) $$','interpreter','latex');grid;
%     legend('0.015','0.02','0.03','0.04','0.05');
    hold on; plot(AB,get(gca,'ylim'),'--y')
%     print -f1 -r300 -depsc2 ../Figures/rga_figure_2.eps
% print -f1 -r300 -depsc2 ../Figures/rga_figure.eps
legend 'rga' 'AB=0.5'

