%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        GLOBAL VARIABLES FILE
%
% Defines global variables necessary to define cost functional JN and
% gradient dJN. We use a MATLAB optimization procedure (lsqnonneg or
% fmincon) which takes JN and dJN as arguments. The optimization 
% procedure requires that the only arguments of JN and dJN are the
% arguments to be optimized. Therefore we define the known parameters of JN
% and dJN as global variables.

load('data\preprocessed.mat','n','tau','nSPLHR','N','T','PAD')


global N n tau T nSPLHR

global P
% P+1 splines/episode * 1 ep/n steps * 1step/tau hr = (P+1)/(tau*n) splines per hour
% We'll use nSPLHR splines per hour, giving
P=round(n*nSPLHR*tau)-1;

N=8; % state and evolution operator discretization index  
%T = n*tau;

global M_state L R K_state
    % LINEAR SPLINE MATRIX
    M_state = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    K_state = N * gallery('tridiag',-ones(1,N),[1,2*ones(1,N-1),1],-ones(1,N));
    L = diag([zeros(1,N),1]);
    R = diag([1,zeros(1,N)]);
    

global SplinesP_linear derivsP_linear
% (P+1) x n linear spline matrix
%tt=0:1/n:(1-(1/n));


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
global q_init c_init parms_init lambda1 lambda2 lambda3
    q_init = [0.85,0.5];
    c_init = 0*[ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
    %lambda1 = 3e-2; % these values good for real data I guess
    %lambda2 = 2e-3;
    lambda1 = 2e-4;
    lambda2 = 1e-5;
    parms_init = [q_init,c_init];
    lambda3 = 1e-4;
    
% Regularization
global Reg dReg
    Reg = @(q,c) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)' + lambda3*(q(1)+q(2));
    dReg = @(q,c) [lambda3,lambda3, 2*c*(lambda1*(SplinesP_linear*SplinesP_linear')+lambda2*(derivsP_linear*derivsP_linear'))];%,ctou(c)*ctou(c)',ctodu(c)*ctodu(c)'];