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

global M

M = 1;

if M <= 0
    M=0;
    globals
    return
end


SplineHandles = cell(1,M+1);
for k=0:M
    SplineHandles{k+1} = @(x) (x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*x+k+1);
end



T = n*tau;

global M_state L R K_state
    % LINEAR SPLINE MATRIX
    M_state = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    L = diag([zeros(1,N),1]);
    R = diag([1,zeros(1,N)]);
    K_state = zeros(N+1,N+1,M+1);

    
global si_diag si_offdiag
    % si short for "spline integrals"
    % each, when multiplied on the left by q1M, gives the diag, lower,
    % upper of Kq, respectively.
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=integral(SplineHandles{k},max(j-2,0)/N,min(j,N)/N);
            si_offdiag(k,j)=integral(SplineHandles{k},max(j-1,0)/N,j/N);
        end
        si_diag(k,N+1)=integral(SplineHandles{k},(N-1)/N,1);
        %K_state(:,:,k) = N*(diag(-si_offdiag(k,:),-1)+diag(si_diag(k,:))+diag(-si_offdiag(k,:),1));
    end

global SplinesP_linear derivsP_linear
% (P+1) x n linear spline matrix


SplinesP_linear = zeros(P+1,n);
derivsP_linear = zeros(P+1,n);

for k=0:P
    linearfun = @(x)(x>(k-1)/P).*(x<=k/P).*(P*x-(k-1)) + (x>k/P).*(x<(k+1)/P).*(-P*x+k+1);
    for i = 1:n
        SplinesP_linear(k+1,i) = integral(linearfun,(i-1)/n,i/n)*n;
        derivsP_linear(k+1,i) = (T*(k-1)/(tau*P)<=i)*(i<T*k/(tau*P))*P/T +(T*(k)/(tau*P)<=i)*(i<T*(k+1)/(tau*P))*-P/T;
    end
end

global dAN_dq
    dAN_dq = zeros(N+1,N+1,M+2);
    % Last page is q2, all zeros.
    for k=1:M+1
        dAN_dq(:,:,k) = -M_state\K_state(:,:,k);
    end

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
global q_init c_init parms_init
    q_init = ones(1,M+2);
    c_init = 1*[ones(1,fix(P/3)),zeros(1,P+1-fix(P/3))];
    lambda1 = 4e-5;
    lambda2 = 1e-6;
    parms_init = [q_init,c_init];

    
% Regularization
global Reg dReg
    Reg = @(q,c) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)';
    dReg = @(q,c) [zeros(1,2), 2*c*(lambda1*(SplinesP_linear*SplinesP_linear')+lambda2*(derivsP_linear*derivsP_linear'))];%,ctou(c)*ctou(c)',ctodu(c)*ctodu(c)'];