
global NUM_EPISODES nSPLHR N M
NUM_EPISODES = 4;
kkkk = fix(sqrt(NUM_EPISODES));
llll = ceil(NUM_EPISODES/kkkk);


% n = 100-2*PAD;
% nn = 100;
% T = 5;
% maxdilation = 2;
% 
% tau = T/(n-1);
% t = 0:tau:T;
%M = 7;
nSPLHR = 2;
N = 128;
M = 32;
%preprocess
%globals
save('data\preprocessed.mat','n','tau','T','nSPLHR','N','M','PAD','NUM_EPISODES','kkkk','llll')
globals_qq

q = zeros(NUM_EPISODES,34);
% WHEN ALLOWING FOR VARYING DIFFUSIVITY CHANGE SIZE OF q

% QUADRATIC PARAMETER
xx = 0:(1/32):1;
q(:,1:end-1) = 4*repmat(abs((xx-.5).^2),NUM_EPISODES,1)+.1;
q(:,end) = .3;
% % LINEARLY VARYING PARAMETER
% q(:,1) = round(abs(1 + .1*randn(1,NUM_EPISODES)),1);
% q(:,2) = round(abs(.3 + .1*randn(1,NUM_EPISODES)),1);
% q(:,3) = round(abs(.3 + .1*randn(1,NUM_EPISODES)),1);
% % CONSTANT PARAMETER
% q(:,1) = 1;%round(abs(1 + .2*randn(1,NUM_EPISODES)),1);
% q(:,2) = .3;%round(abs(.3 + .1*randn(1,NUM_EPISODES)),1);


y_total = mymodel_qq(q,u_total);
y_total = abs(y_total+0.04*y_total.*randn(size(y_total)));

u_total = u_total(1:NUM_EPISODES,1:145);
y_total = y_total(1:NUM_EPISODES,1:145);

t_out = 0:tau:(tau*(size(u_total,2)-1));
n = length(t_out);
T = t_out(end);

N = 8;
save('data\preprocessed.mat','u_total','y_total','t_out','n','N','T','tau','nSPLHR','PAD','NUM_EPISODES','kkkk','llll','q','M')

plot(t_out,u_total')
hold on
plot(t_out,y_total')