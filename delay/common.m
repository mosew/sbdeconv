%% Constant matrices for delay system


global N T n
N = 8;
T = 10;
n = 4*N+1;

lambda1 = 1e-2;
lambda2 = 1e-2;
lambda3 = 1e-2;

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
        dUSplLin(i,k+1) = (T*(k-1)/(tau*N)<=i)*(i<T*k/(tau*N))*N/T +(T*(k)/(tau*N)<=i)*(i<T*(k+1)/(tau*N))*-N/T;
    end
end

ctou = @(c) USplLin*c;
ctodu = @(c) dUSplLin*c;

global Reg dReg
Reg = @(q,c) lambda1*ctou(c)*ctou(c)' + lambda2*ctodu(c)*ctodu(c)' + lambda3*(q(1)+q(2)+q(3)+q(4));
dReg = @(q,c) [lambda3,lambda3,lambda3,lambda3, 2*c*(lambda1*(USplLin*USplLin')+lambda2*(dUSplLin*dUSplLin'))];
