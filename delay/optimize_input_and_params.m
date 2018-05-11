%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Main optimization routine
%

function [a0_star,a1_star,b1_star,ell_star,u_star] = optimize_input_and_params(training_u,training_y,test_y)

global n ctou
global c_init parms_init

% All begin at t=0.
m = size(training_y,1);

assert(all(size(training_y)==[m n]));
assert(all(size(training_u)==[m n]));
assert(all(size(test_y)==[1 n]));


% Options for constrained minimization
Aeq = [];
Beq = [];
Aineq = [];
Bineq = [];
LB = [-10,-10,0,0,-eps*ones(size(c_init))];
UB = [0,0,10,10,1000*ones(size(c_init))];
NONLCON = [];
OPTIONS = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,...
    'ScaleProblem',true,...
    'Diagnostics','off',...
    'StepTolerance',1e-18,...
    'CheckGradients',false);

% Constrained minimization        
JdJ = @(parms) JdJ0(parms,training_u,[training_y;test_y]);
[parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon(JdJ,parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,NONLCON,OPTIONS);

        
% Collect optimal parameters and deconvolved signal
a0_star = parms_star(1);
a1_star = parms_star(2);
b1_star = parms_star(3);
ell_star = parms_star(4);
u_star  = ctou(parms_star(5:end));

end