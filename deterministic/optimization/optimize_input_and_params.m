%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Main optimization routine
%

function [q1_star,q2_star,u_star,lambda1_star,lambda2_star] = optimize_input_and_params(training,test,P)

global SplinesP_linear PAD n

load('data\preprocessed.mat','y_total','u_total','n','PAD')

global training_u Y
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

% % Set initial parameters
% Moved to globals
% q_init = ones(1,2);
% c_init = [ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
% lambda1_init = 5e-2;
% lambda2_init = 5e-2;
% parms_init = [q_init,c_init,lambda1_init,lambda2_init];
global q_init c_init lambda1_init lambda2_init parms_init

% Options for constrained minimization
% Fix first and last linear splines to be 0.
fixzeros = fix(PAD/n*(P+1));


if size(q_init,2)>2
    global M
    Aeq = [zeros(1,M+1),0,1,zeros(1,P-fixzeros),ones(1,fixzeros),0,0];
else
    M=0;
    Aeq = [0,0,1,zeros(1,P-fixzeros),ones(1,fixzeros),0,0];
end

Beq = 0;
Aineq = [];
Bineq = [];
LB = [eps*ones(size(q_init)),-eps*ones(size(c_init)),0.9*lambda1_init,0.9*lambda2_init];
UB = [inf*ones(size(q_init)),inf*ones(size(c_init)),1.1*lambda1_init,1.1*lambda2_init];
NONLCON = [];
OPTIONS = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,...
    'ScaleProblem','none',...
    'MaxIterations',2000,...
    'StepTolerance',eps,...
    'Diagnostics','off',...
    'CheckGradients',false);

% Constrained minimization        
if length(q_init)>2
    [parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon(@J_gradJ_qq,parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
else
    [parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon(@J_gradJ,parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
end

        
% Collect optimal parameters and deconvolved signal
q1_star = parms_star(1:M+1);
q2_star = parms_star(M+2);
u_star  = parms_star(M+3:(end-2))*SplinesP_linear;
u_star = u_star(1:n-PAD);
lambda1_star = parms_star(end-1);
lambda2_star = parms_star(end);

end