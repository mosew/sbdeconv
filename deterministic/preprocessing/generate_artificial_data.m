
global PAD NUM_EPISODES nSPLHR N
PAD = 9;
NUM_EPISODES = 9;
kkkk = fix(sqrt(NUM_EPISODES));
llll = ceil(NUM_EPISODES/kkkk);


n = 100-2*PAD;
nn = 100;
T = 4;
maxdilation = 2;
tau = T/(n-1);
t = 0:tau:T;
nSPLHR = 3;
N = 128;

save('data\preprocessed.mat','n','tau','T','nSPLHR','N','PAD','NUM_EPISODES','kkkk','llll')
globals

rng(7)
f = @(x) abs((randi(2)-1)*0.7*abs(randn+1)*betapdf(x/4,8+randi(4),5+randi(4)) + 0.8*abs((randn+1))*betapdf(x/4,5,9+randi(4)));
u_total = zeros(NUM_EPISODES,nn*maxdilation);

% WHEN ALLOWING FOR VARYING DIFFUSIVITY CHANGE SIZE OF q
q = zeros(NUM_EPISODES,3);
%q(:,1) = abs(1 + trandn(-ones(1,NUM_EPISODES),inf*ones(1,NUM_EPISODES)));
q(:,1) = abs(.5 + .1*randn(1,NUM_EPISODES))+.3;
q(:,2) = abs(.3 + .1*randn(1,NUM_EPISODES))+.3;
%q(:,2) = abs(.3 +.3*trandn(-ones(1,NUM_EPISODES),inf*ones(1,NUM_EPISODES)));
q(:,3) = abs(.3 + .05*randn(1,NUM_EPISODES))+.1;

for i = 1:NUM_EPISODES
    u_total(i,1:nn) = [zeros(1,PAD),f(t),zeros(1,PAD)];
end
y_total = mymodel(q,u_total);
y_total = abs(y_total+0.04*y_total.*randn(size(y_total)));

t_out = 0:tau:(tau*(size(u_total,2)-1));
n = length(t_out);

save('data\preprocessed.mat','u_total','y_total','t_out','n','N','T','tau','nSPLHR','PAD','NUM_EPISODES','kkkk','llll','q')

plot(t_out,u_total')
hold on
plot(t_out,y_total')