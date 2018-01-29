%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Main optimization routine
%

function [q1_star,q2_star,u_star,lambda1_star,lambda2_star] = optimize_input_and_params(training,test,P)

global SplinesP_linear PAD

load('data\preprocessed.mat','y_total','u_total','n','PAD')

global m training_u Y
% All begin at t=0.
training_y = y_total(training,:);
training_u = u_total(training,:);
test_u = u_total(test,:);
test_y = y_total(test,:);
m = size(training_y,1);
Y=[training_y;test_y];


assert(all(size(training_y)==[m n]));
assert(all(size(training_u)==[m n]));
assert(all(size(test_y)==[1 n]));
assert(all(size(test_u)==[1 n]));

% Set initial parameters
q_init = ones(1,2);
c_init = [ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
lambda1_init = 5e-2;
lambda2_init = 5e-2;
parms_init = [q_init,c_init,lambda1_init,lambda2_init];

% Options for constrained minimization
% Fix first and last linear splines to be 0.
fixzeros = fix(n/(P*PAD));
Aeq = [0,0,1,zeros(1,P-fixzeros),ones(fixzeros),0,0];
Beq = 0;
Aineq = [];
Bineq = [];
LB = [eps*ones(size(q_init)),-eps*ones(size(c_init)),1e-9,1e-9];
UB = [inf*ones(size(q_init)),inf*ones(size(c_init)),10,10];
NONLCON = [];
OPTIONS = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,...
    'ScaleProblem','none',...
    'MaxIterations',1000,...
    'StepTolerance',eps,...
    'Diagnostics','on',...
    'CheckGradients',false);

% Constrained minimization        
[parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon(@JN_and_dJN_uspline,parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,NONLCON,OPTIONS);

        
% Collect optimal parameters and deconvolved signal
q1_star = parms_star(1);
q2_star = parms_star(2);
u_star  = parms_star(3:(end-2))*SplinesP_linear;
u_star = u_star(1:n-PAD);
lambda1_star = parms_star(end-1);
lambda2_star = parms_star(end);

end