
function [Bc11,Bc22]=calcula_bounds_nuevo(w,P_som,D_mat,Bs,mB,Xs,phs)

nplants=size(P_som,3);
P110=1/P_som(1,1,1);
% P220=1/P_som(2,2,1);

b11=mB(1,1)*Bs(1,1);
% b21=mB(2,1)*Bs(2,1);
x11=Xs(1,1);
% x21=Xs(2,1);

bc11=b11*x11;
% bc21=b21*x21; 

QD11 = minreal(-P_som(1,1)*D_mat(1,1)-P_som(1,2)*D_mat(2,1));
% QD21 = minreal(-P_som(2,1)*D_mat(1,1)-P_som(2,2)*D_mat(2,1));

aux11=freqresp(P_som(1,1),w);
aux12=freqresp(P_som(1,2),w);
% aux21=freqresp(P_som(2,1),w);
% aux22=freqresp(P_som(2,2),w);

p11_aux=reshape(aux11,size(aux11,3),size(aux11,4)).';
p12_aux=reshape(aux12,size(aux12,3),size(aux12,4)).';
% p21_aux=reshape(aux21,size(aux21,3),size(aux21,4)).';  
% p22_aux=reshape(aux22,size(aux22,3),size(aux22,4)).';
         
b11_aux=repmat(b11,nplants,1);
% b21_aux=repmat(b21,nplants,1); 

d11e_p11=abs(p12_aux).*b21_aux./p11_aux;
% d21e_p22=abs(p21_aux).*b11_aux./p22_aux;

B11_coup=genbnds(10,w,bc11,d11e_p11,0,1,1/P_som(1,1),P110,phs);
% B21_coup=genbnds(10,w,bc21,d21e_p22,0,1,1/P_som(2,2),P220,phs);
B11_distff=sisobnd14c(w,tf(0,1),b11-bc11,1/P_som(1,1),QD11,1,phs);
% B21_distff=sisobnd14c(w,tf(0,1),b21-bc21,1/P_som(2,2),QD21,1,phs);
B11_distff_ad=adaptation(B11_distff,12);
% B21_distff_ad=adaptation(B21_distff,12);

Bc11=grpbnds(B11_coup,B11_distff_ad);
Bc22=[]; %grpbnds(B21_coup,B21_distff_ad);


function plotbndsc(bdb,ptype,phase)


pos=[1 1 1 1];
[bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,ptype,phase,pos);


[rbdb,cbdb]=size(bdb);
lwbs2=length(wbs2);

if lwbs2 ~= cbdb,
    
    a=gca;
    apos=get(a,'pos');
    set(a,'box','on','xgrid','on','ygrid','on',...
        'gridlinestyle',':','drawmode','fast',...
        'nextplot','add','xlim',axs(1:2),'ylim',axs(3:4));
    
    bnd=qplotbd(phase,bdb,coora,coorb,axs);
    
else
    disp('Unplottable bounds.');
end


% Cargar plantas
load plantas

P = sys_red(:,1:2,:);
D = sys_red(:,3,:);
P_som = inv(handles.P);

nplants = size(handles.P,3);
w = [0.01 0.1 1 3 5 10 30 100];
phs = 0:-0.1:-360;

% Especificaciones 
s = tf('s');
b11 = -3600*(s/2.5+1)/(s/1+1)/(s/1+1);
b21 = 4000*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);
b11s(1,:) = squeeze(abs(freqresp(b11,handles.w)));
b21s(1,:) = squeeze(abs(freqresp(b21,handles.w)));


% A partir de aquí se calcula para la frecuencia en concreto
x11_ = 0.3;
i=1;
[Bc11,Bc22] = calcula_bounds_nuevo(w(i),P_som,D,Bs,[1;1],[x11_;0],phs); 
figure(1); plotbndsc(Bc11,[],phs);





for i=1:length(w)
    Bs=[b11s(i); b21s(i)];
    [Bc11,Bc22] = calcula_bounds_nuevo(w(i),P_som,D,Bs,[1;1],[x11_;0],phs); 

    figure(1); plotbndsc(Bc11,[],phs);

