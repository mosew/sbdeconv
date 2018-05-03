%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        GLOBAL VARIABLES FILE
%
% Defines global variables necessary to define cost functional JN and
% gradient dJN. We use a MATLAB optimization procedure (lsqnonneg or
% fmincon) which takes JN and dJN as arguments. The optimization 
% procedure requires that the only arguments of JN and dJN are the
% arguments to be optimized. Therefore we define the known parameters of JN
% and dJN as global variables.

load('data\preprocessed.mat','n','tau','nSPLHR','N','T','PAD','M')


global N n tau T nSPLHR

global P
% P+1 splines/episode * 1 ep/n steps * 1step/tau hr = (P+1)/(tau*n) splines per hour
% We'll use nSPLHR splines per hour, giving
P=round(n*nSPLHR*tau)-1;

SplineHandles = cell(1,M+1);
for k=0:M
    SplineHandles{k+1} = @(x) (k>=0).*(x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*(x-(k+1))).*(k<=M);
end

global M_state L R
    % LINEAR SPLINE MATRIX
    M_state = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    L = diag([zeros(1,N),1]);
    R = diag([1,zeros(1,N)]);
    %K_state = zeros(N+1,N+1,M+1);

    
global si_diag si_offdiag
    % si short for "spline integrals"
    % each, when multiplied on the left by q1M, gives the diag, lower,
    % upper of Kq, respectively.
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=N^2*integral(SplineHandles{k},max(j-2,0)/N,min(j,N)/N);
            si_offdiag(k,j)=N^2*integral(SplineHandles{k},max(j-1,0)/N,j/N);
        end
        si_diag(k,N+1)=N^2*integral(SplineHandles{k},(N-1)/N,1);
        %K_state(:,:,k) = N*(diag(-si_offdiag(k,:),-1)+diag(si_diag(k,:))+diag(-si_offdiag(k,:),1));
    end

global SplinesP_linear derivsP_linear
% (P+1) x n linear spline matrix. Multiplication on the left by c gives
% AVE approximation for linear spline approx. represented by c
SplinesP_linear = zeros(P+1,n);
derivsP_linear = zeros(P+1,n);

for k=0:P
    linearfun = @(x)(k>=0).*(x>(k-1)/P).*(x<=k/P).*(P*x-(k-1)) + (x>k/P).*(x<(k+1)/P).*(-P*x+k+1).*(k<=P);
    for i = 1:n
        SplinesP_linear(k+1,i) = integral(linearfun,(i-1)/n,i/n)*n;
        derivsP_linear(k+1,i) = (T*(k-1)/(tau*P)<=i)*(i<T*k/(tau*P))*P/T +(T*(k)/(tau*P)<=i)*(i<T*(k+1)/(tau*P))*-P/T;
    end
end

M_u = 1/6*(T/P)*gallery('tridiag',ones(1,P),[2,4*ones(1,P-1),2],ones(1,P));
K_u = (P/T)*(diag(-ones(1,P),-1)+diag([1,2*ones(1,P-1),1])+diag(-ones(1,P),1));

global dA_dq
    dA_dq = build_dAN_dqM();

global dB_dq
    dB_dq = zeros(N+1,1,M+2);
    dB_dq(1,:,M+2) = 1;
    dB_dq(:,:,M+2) = M_state\dB_dq(:,:,M+2);
    
global Chat
    Chat = [zeros(1,N),1];

global ctou ctodu
ctou = @(c) c*SplinesP_linear;
ctodu = @(c) c*derivsP_linear;

% Set initial parameters
global q_init c_init parms_init
    %s = 10^((log2(N)-1)/4);
    %q_init = [2*s,1*s,.3];
    q_init = [ones(1,M+1),.3];
    if numel(q_init)>2
        assert(numel(q_init)==M+2);
    end
    c_init = 1*[ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
%     lambda1 = 2e-2;
%     lambda2 = 4e-3;
    lambda1 = 2e-4;
    lambda2 = 1e-5;
    parms_init = [q_init,c_init];

% Regularization
global Reg dReg
    %Reg = @(q,c) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)';
    %dReg = @(q,c) [zeros(1,M+2), 2*c*(lambda1*(SplinesP_linear*SplinesP_linear')+lambda2*(derivsP_linear*derivsP_linear'))];

    %lambda3 = 3e-3;
    lambda3 = 1e-4;

    M_q = 1/(6*M) * gallery('tridiag',ones(1,M),[2,4*ones(1,M-1),2],ones(1,M)); % [<psi^N_i,psi^N_j>]_{ij}
    K_q = M * gallery('tridiag',-ones(1,M),[1,2*ones(1,M-1),1],-ones(1,M));
    Reg = @(q,c) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)' + lambda3*(q(1:end-1)*M_q)*(q(1:end-1)*M_q)';
    dReg = @(q,c) [2*q(1:end-1)*lambda3*(M_q*M_q'),0, 2*c*(lambda1*(SplinesP_linear*SplinesP_linear')+lambda2*(derivsP_linear*derivsP_linear'))];