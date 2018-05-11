common

training = [1,2,4];
test = 3;
[a0_star,a1_star,ell_star,u_star] = optimize_input_and_params(us(training,:),ys(training,:),ys(test,:))

figure
plot(t,us(test,:),'k+');
hold on
plot(t,u_star);

figure
plot(t,ys(test,:),'r+')
hold on
plot(t,delaymodel(a0_star,a1_star,ell_star,u_star,t))