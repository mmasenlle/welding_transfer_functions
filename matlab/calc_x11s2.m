load plantas

P = sys_red(:,1:2,:);
D = sys_red(:,3,:);
P_som = inv(P);

nplants = size(P,3);
w = [0.01 0.1 1 3 5 10 30 100];
phs = 0:-0.1:-360;

% Especificaciones 
s = tf('s');
b11 = -3600*(s/2.5+1)/(s/1+1)/(s/1+1);
b21 = 4000*(s^2/539.4126^2+2*s*1/539.4126+1)/(s/1.0007+1)/(s/1.0007+1);
b11s(1,:) = squeeze(abs(freqresp(b11,w)));
b21s(1,:) = squeeze(abs(freqresp(b21,w)));


% A partir de aqu� se calcula para la frecuencia en concreto
x11_ = 0.35;
i=1;
[Bc11,Bc22] = calcula_bounds(w(i),P_som,D,[b11s(i);b21s(i)],[1;1],[x11_;0],phs);
norm(Bc11(:,2)-Bc11(:,1))
figure(1); plotbndsc(Bc11,[],phs);

x11w=zeros(length(w));
for i=1:length(w)
    Bs=[b11s(i); b21s(i)];
    x11_ = 0.0;
    step_ = 0.1;
    last_x11 = x11_;
    best_v = 1e20;
    for k=1:30
        [Bc11,Bc22] = calcula_bounds(w(i),P_som,D,[b11s(i);b21s(i)],[1;1],[x11_;0],phs);
        v=norm(Bc11(:,2)-Bc11(:,1));
        if v > best_v
            if step_ > 0.0011
                step_ = step_ * 0.1;
                x11_ = last_x11;
            else
                x11w(i) = last_x11;
                break;
            end
        else
            best_v = v;
            last_x11 = x11_; 
        end
        x11_ = x11_ + step_;
    end
end
w
x11w
%figure(1); plotbndsc(Bc11,[],phs);

% 0.47 0.48 *0.21* 0.22 0.23 0.11 0.07 
% 0.1000    0.1100    0.3000    0.1010    0.3000    0.2010         0         0
% 0.1000    0.1100    0.1920    0.3000    0.2000    0.1000    0.1000    0.0140
% 0.5300    0.5500    0.2100         0         0         0         0         0
