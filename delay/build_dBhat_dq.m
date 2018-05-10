function dBhat_dq = build_dBhat_dq(A,dA_dq,Ahat,dAhat_dq,B)
    global N
    dBhat_dq = zeros(2*N+1,1,4);
    
    dBhat_dq(:,1,1) = -A\(dA_dq(:,:,1)/A*(Ahat - eye(N+1)) - dAhat_dq(:,:,1))*B;
    dBhat_dq(:,1,2) = -A\(dA_dq(:,:,2)/A*(Ahat - eye(N+1)) - dAhat_dq(:,:,2))*B;
    dBhat_dq(:,1,4) = -A\(dA_dq(:,:,4)/A*(Ahat - eye(N+1)) - dAhat_dq(:,:,4))*B;
end