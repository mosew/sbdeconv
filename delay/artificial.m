%% Artificial model for delay system
common
ell = 1;
a0 = -0.01;
a1 = -0.3;
b1 = 0.1;

global T n h
t = 0:(T/(n-1)):T;
h = T/(n-1);

k=build_KN(a0,a1,b1);
assert(all(eig(k'*k)~=0));

u1 = @(s) exp(-(s-2).^2);
u2 = @(s) exp(-(s-1.8).^2/.5)+0.4*exp(-(s-4).^2);
u3 = @(s) 0.1*s;
u4 = @(s) ones(size(s));

u1 = u1(t);
u2 = u2(t);
u3 = u3(t);
u4 = u4(t);

y1 = delaymodel(a0,a1,b1,ell,u1,t);
y2 = delaymodel(a0,a1,b1,ell,u2,t);
y3 = delaymodel(a0,a1,b1,ell,u3,t);
y4 = delaymodel(a0,a1,b1,ell,u4,t);

us = [u1;u2;u3;u4];
ys = [y1;y2;y3;y4];

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