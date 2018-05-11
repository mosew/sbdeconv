function [JN,dJN] = JdJ0(parmsu,tru,trY)
    % INPUTS
    % parmsu: 1 x 3+N+1 matrix of a0,a1,ell,linear spline coefficients for unknown u
    % tru: m x n matrix of linear spline coefficients for training u
    % trY: (m+1) x n matrix of sampled outputs

    global ctou N n h
    global Reg dReg USplLin

    [m,~]=size(tru);
    
    
    % Process inputs
    a0 = parmsu(1);
    a1 = parmsu(2);
    ell = parmsu(3);
    u = parmsu(4:end);
    % CAREFUL TO SCALE DATA IF NECESSARY
    
    MN = build_MN(ell);
    
    total_u = [tru;ctou(u)'];

    % DEFINE SYSTEM OPERATORS
    KN = build_KN(a0,a1);
    AN = MN\KN;
    BN = [zeros(1,N+1),2*N/ell,zeros(1,N-1)]';
    CN = [1,zeros(1,2*N)]*MN;

    
    %sg = @(s) expm(s*AN);
    %sgh = sg(h);

    dAN_dq  = zeros(2*N+1,2*N+1,4);
    dAN_dq(:,:,1) = MN\build_KN(1,a1);
    dAN_dq(:,:,2) = MN\build_KN(a0,1);
    dAN_dq(:,:,3) = -MN\(build_MN(1)*AN);
    
    dAhat_dq = zeros(2*N+1,2*N+1,3);

    for k = 1:3
        AdAExp = expm([h*AN, h*dAN_dq(:,:,k);zeros(2*N+1),h*AN]);
        dAhat_dq(:,:,k) = AdAExp(1:(2*N+1),(2*N+2):end);
    end
    Ahat = AdAExp(1:(2*N+1),1:(2*N+1));
    Bhat = (Ahat - eye(2*N+1))*(AN\BN);
    dBhat_dq = build_dBhat_dq(AN,dAN_dq,Ahat,dAhat_dq,BN,ell);

    X = zeros(2*N+1,n,m+1);
    
    for i = 1:m+1
        X(:,1,i) = 0;% assumed u(i,1)=0, and initial state=0;
        for j = 2:n
            X(:,j,i) = Ahat*X(:,j-1,i) + Bhat*total_u(i,j-1);
        end
    end

    % INITIALIZE ETA
    eta = zeros(2*N+1,n,m+1);
    % goes from t=h to t=h*n. At t=0, everything is 0...
    for i = 1:m+1
        % episodes in Y are assumed to start at 0.
        eta(:,n,i) = (CN*X(:,n,i)-trY(i,n))*CN';
    end


    % COMPUTE GRADIENT CONTRIBUTIONS
    dJN = zeros(1,3+N+1); % one for each component of (q,c)
    JN=0;

    % Adjoint method.
    for j=(n-1):-1:1
        for i=1:m+1
            eta(:,j,i) = Ahat' * eta(:,j+1,i) + (CN*X(:,j,i)-trY(i,j))*CN';
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
                for k = 1:3
                    dJN(k) = dJN(k) + sc*(eta(:,j,i)' * (dAhat_dq(:,:,k)*X(:,j-1,i) + dBhat_dq(:,:,k)*total_u(i,j-1)));
                end
            end
            JN = JN + sc*(CN*X(:,j,i)-trY(i,j))^2;
        end

        for r=0:N
            dJN(r+4) = dJN(r+4) + sc*(eta(:,j,m+1)'*Bhat*USplLin(j,r+1));
        end
    end

    JN = JN/m + Reg([a0,a1,ell],u);
    dJN = dJN/m + dReg([a0,a1,ell],u);
    
end