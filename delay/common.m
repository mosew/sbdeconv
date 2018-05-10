%% Constant matrices for delay system


global N T n h
N = 8;
T = 10;
n = 3*N+1;
h = T/(n-1);

lambda1 = 1e-8;
lambda2 = 1e-8;
lambda3 = 1e-8;

global q_init c_init parms_init
q_init = [0.1,-0.1,0.1,0.8];
c_init = zeros(1,N+1);
parms_init = [q_init,c_init];

global scrKN scrKNm scrMNdell scrMNmdell ctou USplLin dUSplLin
scrMNdell = 1/(6*N)*gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N));
scrMNmdell = scrMNdell(2:end,2:end);
scrKN = 0.5*gallery('tridiag',ones(1,N),[1,zeros(1,N-1),-1],-ones(1,N));
scrKNm = scrKN(2:end,2:end);

USplLin = zeros(n,N+1);
dUSplLin = zeros(n,N+1);
for k=0:N
    linearfun = @(x)(x>(k-1)/N).*(x<=k/N).*(N*x-(k-1)) + (x>k/N).*(x<(k+1)/N).*(-N*x+k+1);
    for i = 1:n
        USplLin(i,k+1) = integral(linearfun,(i-1)/n,i/n)*n;
        dUSplLin(i,k+1) = (T*(k-1)/(h*N)<=i)*(i<T*k/(h*N))*N/T +(T*(k)/(h*N)<=i)*(i<T*(k+1)/(h*N))*-N/T;
    end
end

ctou = @(c) USplLin*c';
ctodu = @(c) dUSplLin*c';

global Reg dReg
Reg = @(q,c) lambda1*ctou(c)'*ctou(c) + lambda2*ctodu(c)'*ctodu(c) + lambda3*(q(1)+q(2)+q(3)+q(4));
dReg = @(q,c) [lambda3,lambda3,lambda3,lambda3, 2*c*(lambda1*(USplLin'*USplLin)+lambda2*(dUSplLin'*dUSplLin))];
