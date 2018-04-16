function dBhat_dq = build_dBhat_dq(A,Ahat,dAhat_dq,B)
    global N dA_dq dB_dq
    dBhat_dq = zeros(N+1,1,2);
    
    dBhat_dq(:,1,1) = -A\(dA_dq(:,:,1)/A*(Ahat - eye(N+1)) - dAhat_dq(:,:,1))*B;
    dBhat_dq(:,1,2) =A\(Ahat-eye(N+1))*dB_dq(:,:,2);
end