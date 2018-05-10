%% Artificial model for delay system

ell = 0.6;
a0 = -.2;
a1 = -.5;
b1 = .3;

global T h
t = 0:h:T;

u1 = @(s) exp(-(s-2).^2);
u2 = @(s) exp(-(s-1.8).^2/8)+0.4*exp(-(s-3).^2/4);
u3 = @(s) 0.1*s;
u4 = @(s) ones(size(s));

y1 = delaymodel(a0,a1,b1,ell,u1,t);
y2 = delaymodel(a0,a1,b1,ell,u2,t);
y3 = delaymodel(a0,a1,b1,ell,u3,t);
y4 = delaymodel(a0,a1,b1,ell,u4,t);

us = [u1(t);u2(t);u3(t);u4(t)];
ys = [y1;y2;y3;y4];

subplot(2,2,1)
plot(t,us(1,:),t,y1)
legend('u','y')

subplot(2,2,2)
plot(t,u2(t),t,y2)
legend('u','y')

subplot(2,2,3)
plot(t,u3(t),t,y3)
legend('u','y')

subplot(2,2,4)
plot(t,u4(t),t,y4)
legend('u','y')

training = 1:3;
test = 4;
[a0_star,a1_star,b1_star,ell_star,u_star] = optimize_input_and_params(us(training,:),ys(training,:),ys(test,:))