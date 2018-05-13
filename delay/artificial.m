%% Artificial model for delay system
common
ell = 2;
a0 = -0.5;

global T n h
t = 0:(T/(n-1)):T;
h = T/(n-1);

k=build_KN(a0);
assert(all(eig(k'*k)~=0));

u1 = @(s) exp(-(s-5).^2/5);
u2 = @(s) exp(-(s-3).^2/2)+0.4*exp(-(s-10).^2/4);
u3 = @(s) 0.1*s;
u4 = @(s) ones(size(s));

u1 = u1(t);
u2 = u2(t);
u3 = u3(t);
u4 = u4(t);

y1 = delaymodel(a0,ell,u1,t);
y2 = delaymodel(a0,ell,u2,t);
y3 = delaymodel(a0,ell,u3,t);
y4 = delaymodel(a0,ell,u4,t);

us = [u1;u2;u3;u4]/10;
ys = [y1;y2;y3;y4]/10;
save('artificial.mat','us','ys','N','T')

subplot(2,2,1)
plot(t,us(1,:),t,y1)
legend('u','y')

subplot(2,2,2)
plot(t,u2,t,y2)
legend('u','y')

subplot(2,2,3)
plot(t,u3,t,y3)
legend('u','y')

subplot(2,2,4)
plot(t,u4,t,y4)
legend('u','y')