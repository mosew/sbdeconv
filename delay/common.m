%% Constant matrices for delay system


global N T n h
N = 16;
T = 30;
n = 4*N+1;
h = T/(n-1);

global q_init c_init parms_init
q_init = [-0.06,-0.2,0.04,0.7];
c_init = eps*ones(1,N+1);
parms_init = [q_init,c_init];

lambda1 = 0;%;1;
lambda2 = 0;%0.002;
lambda3 = 1e-10;

global scrKN scrKNm scrMNdell scrMNmdell ctou USplLin dUSplLin
scrMNdell = 1/(6*N)*gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N));
scrMNmdell = scrMNdell(2:end,2:end);
scrKN = 0.5*gallery('tridiag',ones(1,N),[1,zeros(1,N-1),-1],-ones(1,N));
scrKNm = scrKN(2:end,2:end);

USplLin = zeros(n,N+1);
dUSplLin = zeros(n,N+1);
for k=0:N
    linearfun = @(x)(x>(T/N*(k-1))).*(x<=(T/N*k)).*(N/T*(x-(k-1)*T/N)) + (x>(T/N*k)).*(x<(T/N*(k+1))).*(-N/T*(x-(k+1)*T/N));
    for i = 1:n
        USplLin(i,k+1) = linearfun(h*(i-1));
        dUSplLin(i,k+1) = ((h*(i-1))>T/N*(k-1)).*((h*(i-1))<=T/N*k).*N/T-((h*i)>T/N*(k-1)).*((h*i)<=T/N*k).*N/T;
    end
end

ctou = @(c) USplLin*c';
ctodu = @(c) dUSplLin*c';

global Reg dReg
Reg = @(q,c) lambda1*ctou(c)'*ctou(c) + lambda2*ctodu(c)'*ctodu(c) + lambda3*(q(1)+q(2)+q(3)+q(4));
dReg = @(q,c) [lambda3,lambda3,lambda3,lambda3, 2*c*(lambda1*(USplLin'*USplLin)+lambda2*(dUSplLin'*dUSplLin))];
