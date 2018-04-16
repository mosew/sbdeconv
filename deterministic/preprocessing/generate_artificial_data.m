
global PAD NUM_EPISODES nSPLHR
PAD = 9;
NUM_EPISODES = 9;
kkkk = fix(sqrt(NUM_EPISODES));
llll = ceil(NUM_EPISODES/kkkk);

global n tau
n = 100-2*PAD;
T = 10;
tau = T/(n-1);
t = 0:tau:T;
nSPLHR = 2;

save('data\preprocessed.mat','n','tau','nSPLHR','PAD','NUM_EPISODES','kkkk','llll')
globals

rng(6)
f = @(x) abs(0.3*betapdf(x,12+randi(4)-2,7+randi(4)-2) + 0.6*betapdf(x,5,11+randi(4)-2) + 0.08*randn(1,n));
u_total = zeros(NUM_EPISODES,n+2*PAD);

% WHEN ALLOWING FOR VARYING DIFFUSIVITY CHANGE SIZE OF q
q = zeros(NUM_EPISODES,2);
q(:,1) = abs(1 + .1*randn(NUM_EPISODES,1));
q(:,2) = abs(10 + 2*randn(NUM_EPISODES,1));

for i = 1:NUM_EPISODES
    u_total(i,:) = [zeros(1,PAD),f(t),zeros(1,PAD)];
end
y_total = mymodel(q,u_total);
n = 100;

% for i = 1:NUM_EPISODES
%     subplot(kkkk,llll,i)
%     plot(u_total(i,:))
%     hold on
%     plot(y_total(i,:))
% end


save('data\preprocessed.mat','u_total','y_total','n','tau','nSPLHR','PAD','NUM_EPISODES','kkkk','llll')
