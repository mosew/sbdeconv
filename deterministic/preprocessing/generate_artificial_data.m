
global PAD NUM_EPISODES nSPLHR N M
PAD = 9;
NUM_EPISODES = 4;
kkkk = fix(sqrt(NUM_EPISODES));
llll = ceil(NUM_EPISODES/kkkk);


n = 100-2*PAD;
nn = 100;
T = 4;
maxdilation = 3;

tau = T/(n-1);
t = 0:tau:T;
M = 7;
nSPLHR = 3;
N = 128;

save('data\preprocessed.mat','n','tau','T','nSPLHR','N','PAD','NUM_EPISODES','kkkk','llll')
globals_qq

rng(666)
f = @(x) abs((randi(2)-1)*0.7*abs(randn+1)*betapdf(x/(4*(1+0.3*randn)),8+randi(4),5+randi(4)) + 0.8*abs(randn+1)*betapdf(x/(4*(1+0.3*randn)),5,9+randi(4)));
u_total = zeros(NUM_EPISODES,nn*maxdilation);

% WHEN ALLOWING FOR VARYING DIFFUSIVITY CHANGE SIZE OF q
% QUADRATIC PARAMETER
q = zeros(NUM_EPISODES,M+2);
xx = 0:(1/M):1;
q(:,1:end-1) = 1.*abs((xx-.5).^2+.008*randn(NUM_EPISODES,M+1));
q(:,end) = 1.*abs(.3 + .05*randn(NUM_EPISODES,1))+.1;
% LINEARLY VARYING PARAMETER
%q(:,1) = abs(30 + .1*randn(1,NUM_EPISODES))+.3;
%q(:,2) = abs(10 + .1*randn(1,NUM_EPISODES))+.3;
%q(:,3) = abs(.3 + .05*randn(1,NUM_EPISODES))+.1;
% CONSTANT PARAMETER
%q(:,1) = abs(.7 + .1*randn(1,NUM_EPISODES))+.3;
%q(:,2) = abs(.3 + .05*randn(1,NUM_EPISODES))+.1;

for i = 1:NUM_EPISODES
    u_total(i,1:nn) = [zeros(1,PAD),f(t),zeros(1,PAD)];
end
y_total = mymodel_qq(q,u_total);
y_total = abs(y_total+0.04*y_total.*randn(size(y_total)));

t_out = 0:tau:(tau*(size(u_total,2)-1));
n = length(t_out);
T = t_out(end);

save('data\preprocessed.mat','u_total','y_total','t_out','n','N','T','tau','nSPLHR','PAD','NUM_EPISODES','kkkk','llll','q','M')

plot(t_out,u_total')
hold on
plot(t_out,y_total')