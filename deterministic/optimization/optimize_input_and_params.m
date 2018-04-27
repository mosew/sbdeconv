%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Main optimization routine
%

function [q1_star,q2_star,u_star] = optimize_input_and_params(training,test,P)

global SplinesP_linear nSPLHR n 

load('data\preprocessed.mat','y_total','u_total','n','PAD','nSPLHR')

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

global q_init c_init parms_init

% Options for constrained minimization
% Fix first spline and last hour to be 0.
%fixzeros = fix(PAD/n*(P+1));
fixzeros=nSPLHR;


if size(q_init,2)>2
    M = size(q_init,2)-2;
    Aeq = diag([zeros(1,M+1),0,1,zeros(1,P-fixzeros),ones(1,fixzeros)]);
else
    M=0;
    Aeq = diag([0,0,1,zeros(1,P-fixzeros),ones(1,fixzeros)]);
end

%Aeq=[];
Beq = zeros(size(Aeq,1),1);
Aineq = [];
Bineq = [];
LB = [eps*ones(size(q_init)),-eps*ones(size(c_init))];
UB = [inf*ones(size(q_init)),inf*ones(size(c_init))];
NONLCON = [];
OPTIONS = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,...
    'ScaleProblem','none',...
    'Diagnostics','off',...
    'StepTolerance',1e-16,...
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
u_star  = parms_star((M+3):end)*SplinesP_linear;

end