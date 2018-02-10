function dBhat_dq = build_dBhat_dq_qq(A,Ahat,dAhat_dq,B)
    global N dA_dq dB_dq M
    dBhat_dq = zeros(N+1,1,M+2);
    
    for k = 1:M+1
        dBhat_dq(:,1,k) = -A\(dA_dq(:,:,k)/A*(Ahat - eye(N+1)) - dAhat_dq(:,:,k))*B;
    end
    dBhat_dq(:,1,2) =A\(Ahat-eye(N+1))*dB_dq(:,:,2);
end