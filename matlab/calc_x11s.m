clear all
load plantas

P = sys_red(:,1:2,:);
D = sys_red(:,3,:);
P_som = inv(P);

nplants = size(P,3);
w = [0.01 0.1 1 3 5 10 30 100];
phs = 0:-0.1:-360;

% Especificaciones 
s = tf('s');
% b11 = -3600*(s/2.5+1)/(s/1+1)/(s/1+1);
% b21 = 4000*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);
b11 = -3000*(s/2.5+1)/(s/1+1)/(s/1+1);
b21 = 4200*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);
% b11 = 1311.2001*(s/-0.061672+1)*(s^2/0.55206^2+2*s*0.35452/0.55206+1)*(s^2/0.99684^2+2*s*0.24464/0.99684+1)*(s^2/2.6306^2+2*s*0.29699/2.6306+1)*(s^2/4.4836^2+2*s*0.75119/4.4836+1)/(s^2/0.40632^2+2*s*0.91788/0.40632+1)/(s^2/0.62142^2+2*s*0.53923/0.62142+1)/(s^2/0.99206^2+2*s*0.27172/0.99206+1)/(s^2/1.7454^2+2*s*0.83221/1.7454+1)/(s^2/3.5534^2+2*s*0.60652/3.5534+1);
% b21 = 5993.6235*(s/-56.2953+1)*(s/46.8596+1)*(s^2/0.46234^2+2*s*0.37125/0.46234+1)*(s^2/0.8656^2+2*s*0.27389/0.8656+1)/(s/5.4359+1)/(s/0.38756+1)/(s^2/0.48861^2+2*s*0.68364/0.48861+1)/(s^2/0.81808^2+2*s*0.28174/0.81808+1);
b11s(1,:) = squeeze(abs(freqresp(b11,w)));
b21s(1,:) = squeeze(abs(freqresp(b21,w)));


% % A partir de aqu� se calcula para la frecuencia en concreto
% x11_ = 0.6782;
% i=3;
% Bc11 = calcula_bounds(w(i),P_som,D,[b11s(i);b21s(i)],[1;1],[x11_;0],phs);
% norm(Bc11(:,2)-Bc11(:,1))
% cost_func(Bc11)
% figure(1); plotbndsc(Bc11,[],phs);
% % x11 = [0.4243 0.4503 0.6782 0.4575 0.1976 0.1292 0.0557 0.3777];

% func11 = @(x) cost_func(calcula_bounds(w(i),P_som,D,[b11s(i);b21s(i)],[1;1],[x;0],phs))
% best_x11 = fminbnd(func11,0,1)

% func11(0.5)
% func11(0.999)

x11w=zeros(1,length(w));
for i=1:length(w)
    func11 = @(x) cost_func(calcula_bounds(w(i),P_som,D,[b11s(i);b21s(i)],[1;1],[x;0],phs))
    x11w(i) = fminbnd(func11,0,1);
    x11w
end
w
x11w