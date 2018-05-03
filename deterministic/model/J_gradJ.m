%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Cost functional and its gradient
%
%

function [JN,dJN] = J_gradJ(qc)
 

% x(j+1) = Ax(j) + Bu(j)
% y(j)   = Cx(j)
% x(0)   = 0

% where each x(j) is an N+1-vector representing a state.
% and the operators A,B,C depend on a parameter q = (q1,q2).
% 
%
% Write
%
% u = phi*c
%
% where c=[c(1) .... c(P)]' is a P-vector and
% phi is a rank P, n x P matrix.
%       phi represents the AVE approximations of the linear splines on a
%       P-mesh of the interval [0,n*tau].
% dphi is the AVE approximations of the derivatives of the linear splines.

% Then
% y = L(q,u) = L(q,phi*c) := [sum_{k=0}^{j-1} CA^{j-k-1}B*(phi*c)]_{j=1}^n



% MINIMIZATION PROBLEM:

% We observe paired training data
% (u(i,j),y(i,j)) at times j=1,...,n
% for separate episodes i=1,..,m.

% Given this paired training data and a potentially unseen
% output signal y=[y(1),...,y(n)], determine the inputs u(j).

% J(q(1),q(2),c(1),...c(p),lambda1,lambda2) = 
%     1/m * sum_{i=1}^m sum_{j=1}^n (L(q,u(i,j)) - y(i,j))^2 +
%                       sum_{j=1}^n (L(q, phi*c)  - y0(j))^2
%                + lambda1* ||phi*c||^2 
%                + lambda2* ||dphi*c||^2.


% INPUT:
% qc(1) is q1
% qc(2) is q2
% qc(3:P+3) is c, the linear spline coefficients of the input
% qc(P+4) is lambda1, the regularization coefficient for u
% qc(P+5) is lambda2, the regularization coefficient for u'

% OUTPUT:
% JN is (scalar) value of cost function JN evaluated at qu and
% dJN, the gradient vector of JN, a P+5 - vector

%%%%%

% We have to make a bunch of variables global so that this function, which
% we pass to the optimization toolbox, doesn't take any of the below to be variables to be optimized;
% we could define them here but then they'll be created every
% time J is called.

global N P K_state n M_state L tau dA_dq m
global training_u Y Reg dReg
global Chat SplinesP_linear

% Process inputs
q1 = qc(1);
q2 = qc(2);

% Scale system parameters
c = qc(3:end);
assert(length(c)==P+1);
testu = @(c) c*SplinesP_linear;
total_u = [training_u;testu(c)];

m = size(training_u,1);

% DEFINE SYSTEM OPERATORS
A = -M_state\(L+q1*K_state);
B = M_state\[q2;zeros(N,1)];


mm = [tau*A, tau*dA_dq(:,:,1);zeros(N+1),tau*A];
AdAExp = expm(mm);
Ahat = AdAExp(1:(N+1),1:(N+1));
Bhat = (Ahat - eye(N+1))*(A\B);
dAhat_dq = zeros(N+1,N+1,2);
dAhat_dq(:,:,1) = AdAExp(1:(N+1),(N+2):end);

dBhat_dq = build_dBhat_dq(A,Ahat,dAhat_dq,B);



X = zeros(N+1,n,m+1);

for i = 1:m+1
    X(:,1,i) = 0;%Bhat*total_u(i,1); assumed u(i,1)=0.
    for j = 2:n
        X(:,j,i) = Ahat*X(:,j-1,i) + Bhat*total_u(i,j-1);
    end
end

% INITIALIZE ETA
eta = zeros(N+1,n,m+1);
% goes from t=tau to t=tau*n. At t=0, everything is 0.
for i = 1:m+1
    % episodes in Y are assumed to start at 0.
    eta(:,n,i) = (Chat*X(:,n,i)-Y(i,n))*Chat';
end


% COMPUTE GRADIENT CONTRIBUTIONS
dJN = zeros(1,2+P+1); % one for each component of (q,c)
JN=0;

% Adjoint method.
for j=(n-1):-1:1
    for i=1:m+1
        eta(:,j,i) = Ahat' * eta(:,j+1,i) + (Chat*X(:,j,i)-Y(i,j))*Chat';
    end
end

% Gradient
for j = 1:n
    for i = 1:m+1
        % Gradient of system parameters

        % Scaling constant to decrease importance of system parameters
        % relative to deconvolution

        sc = m*(i==m+1)+(i<=m);
        if j>1
            for k = 1:2
                dJN(k) = dJN(k) + sc*(eta(:,j,i)' * (dAhat_dq(:,:,k)*X(:,j-1,i) + dBhat_dq(:,:,k)*total_u(i,j-1)));
            end
        end
        JN = JN + sc*(Chat*X(:,j,i)-Y(i,j))^2;
    end
    
    for r=0:P
        dJN(r+3) = dJN(r+3) + sc*(eta(:,j,m+1)'*Bhat*SplinesP_linear(r+1,j));
    end
end

JN = JN/m + Reg([q1,q2],c);
dJN = dJN/m + dReg([q1,q2],c);

end