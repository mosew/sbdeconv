%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        GLOBAL VARIABLES FILE
%
% Defines global variables necessary to define cost functional JN and
% gradient dJN. We use a MATLAB optimization procedure (lsqnonneg or
% fmincon) which takes JN and dJN as arguments. The optimization 
% procedure requires that the only arguments of JN and dJN are the
% arguments to be optimized. Therefore we define the known parameters of JN
% and dJN as global variables.

load('data\preprocessed.mat')


global N n tau
N=32; % state and evolution operator discretization index  
tau = 5/60;
T = n*tau;

global M_state L R K_state
    % LINEAR SPLINE MATRIX
    M_state = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    K_state = N * gallery('tridiag',-ones(1,N),[1,2*ones(1,N-1),1],-ones(1,N));
    L = diag([zeros(1,N),1]);
    R = diag([1,zeros(1,N)]);
    
global P nSPLHR
% P+1 splines/episode * 1 ep/n steps * 1step/5min * 60min/1hr = 60(P+1)/(5n) = 12(P+1)/n splines per hour
% We'll use 1 spline per hour, giving P = n/12 - 1
P=round(n*nSPLHR/12)-1;

global SplinesP_linear derivsP_linear
% (P+1) x n linear spline matrix
t=0:1/n:(1-(1/n));


SplinesP_linear = zeros(P+1,n);
derivsP_linear = zeros(P+1,n);

for k=0:P
    linearfun = @(x)(x>(k-1)/P).*(x<=k/P).*(P*x-(k-1)) + (x>k/P).*(x<(k+1)/P).*(-P*x+k+1);
    for i = 1:n
        SplinesP_linear(k+1,i) = integral(linearfun,(i-1)/n,i/n)*n;
        derivsP_linear(k+1,i) = (T*(k-1)/(tau*P)<=i)*(i<T*k/(tau*P))*P/T +(T*(k)/(tau*P)<=i)*(i<T*(k+1)/(tau*P))*-P/T;
    end
end

global dA_dq
    dA_dq = zeros(N+1,N+1,2);
    dA_dq(:,:,1) = -M_state\K_state;

global dB_dq
    dB_dq = zeros(N+1,1,2);
    dB_dq(1,:,2) = 1;
    dB_dq(:,:,2) = M_state\dB_dq(:,:,2);
    
global Chat
    Chat = [zeros(1,N),1];

global ctou ctodu
ctou = @(c) c*SplinesP_linear;
ctodu = @(c) c*derivsP_linear;

% Set initial parameters
global q_init c_init lambda1_init lambda2_init parms_init
    q_init = ones(1,2);
    c_init = [ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
    lambda1_init = 5e-2;
    lambda2_init = 5e-2;
    parms_init = [q_init,c_init,lambda1_init,lambda2_init];

    
% Regularization
global Reg dReg
    Reg = @(q,c,lambda1,lambda2) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)';
    dReg = @(q,c,lambda1,lambda2) [zeros(1,2), 2*c*(lambda1*(SplinesP_linear*SplinesP_linear')+lambda2*(derivsP_linear*derivsP_linear')),ctou(c)*ctou(c)',ctodu(c)*ctodu(c)'];