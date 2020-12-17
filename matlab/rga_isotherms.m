clear all; clc; close all;

load('ss_vars_isotherms')


sys1200 = [ss(A,B,C_1200x,0) ss(A,B,C_1200y,0)
           ss(A,Bv,C_1200x,0) ss(A,Bv,C_1200y,0)];
sys1000 = [ss(A,B,C_1000x,0) ss(A,B,C_1000y,0)
           ss(A,Bv,C_1000x,0) ss(A,Bv,C_1000y,0)];
sys800 = [ss(A,B,C_800x,0) ss(A,B,C_800y,0)
           ss(A,Bv,C_800x,0) ss(A,Bv,C_800y,0)];
sys600 = [ss(A,B,C_600x,0) ss(A,B,C_600y,0)
           ss(A,Bv,C_600x,0) ss(A,Bv,C_600y,0)];


figure
bode(sys1200(1,1),sys1000(1,1),sys800(1,1),sys600(1,1),sys1200(2,2),sys1000(2,2),sys800(2,2),sys600(2,2))
title('Bodes isotherms systems')
legend('1200_1','1000_1','800_1','600_1','1200_2','1000_2','800_2','600_2')
grid


figure
step(sys1200,sys1000,sys800,sys600)
title('Steps isotherms systems')
legend('1200','1000','800','600');
grid

close all
w = logspace(-3,2,100);
Pf1 = freqresp(sys1200,w);
Pf2 = freqresp(sys1000,w);
Pf3 = freqresp(sys800,w);
Pf4 = freqresp(sys600,w);

for i=1:length(w)
    RGAw1(:,:,i)=Pf1(:,:,i).*inv(Pf1(:,:,i)).';
    RGAnum1(:,:,i) = sum(sum(abs(RGAw1(:,:,i) - eye(2))));
    RGAw2(:,:,i)=Pf2(:,:,i).*inv(Pf2(:,:,i)).';
    RGAnum2(:,:,i) = sum(sum(abs(RGAw2(:,:,i) - eye(2))));
    RGAw3(:,:,i)=Pf3(:,:,i).*inv(Pf3(:,:,i)).';
    RGAnum3(:,:,i) = sum(sum(abs(RGAw3(:,:,i) - eye(2))));
    RGAw4(:,:,i)=Pf4(:,:,i).*inv(Pf4(:,:,i)).';
    RGAnum4(:,:,i) = sum(sum(abs(RGAw4(:,:,i) - eye(2))));
end

figure
subplot(2,2,1); semilogx(w,abs(squeeze(RGAw1(1,1,:))),w,abs(squeeze(RGAw2(1,1,:))),w,abs(squeeze(RGAw3(1,1,:))),w,abs(squeeze(RGAw4(1,1,:))))
xlim([1e-3,1e1]);ylabel('RGA');
title('$$ RGAw(P,T_{1}) $$','interpreter','latex');grid;

subplot(2,2,2); semilogx(w,abs(squeeze(RGAw1(1,2,:))),w,abs(squeeze(RGAw2(1,2,:))),w,abs(squeeze(RGAw3(1,2,:))),w,abs(squeeze(RGAw4(1,2,:))))
 xlim([1e-3,1e1]);
title('$$ RGAw(P,T_{2}) $$','interpreter','latex');grid;
legend('1200','1000','800','600');

subplot(2,2,3); semilogx(w,abs(squeeze(RGAw1(2,1,:))),w,abs(squeeze(RGAw2(2,1,:))),w,abs(squeeze(RGAw3(2,1,:))),w,abs(squeeze(RGAw4(2,1,:))))
 xlim([1e-3,1e1]);xlabel('Frecuency (rad/s)');ylabel('RGA')
title('$$ RGAw(V,T_{1}) $$','interpreter','latex');grid;

subplot(2,2,4); semilogx(w,abs(squeeze(RGAw1(2,2,:))),w,abs(squeeze(RGAw2(2,2,:))),w,abs(squeeze(RGAw3(2,2,:))),w,abs(squeeze(RGAw4(2,2,:))))
 xlim([1e-3,1e1]);xlabel('Frecuency (rad/s)')
title('$$ RGAw(V,T_{2}) $$','interpreter','latex');grid;



