function dBhat_dq = build_dBhat_dq(A,dA_dq,Ahat,dAhat_dq,B,ell)
    global N
    dBhat_dq = zeros(2*N+1,1,3);
    
    dBhat_dq(:,1,1) = -A\(dA_dq(:,:,1)/A*(Ahat - eye(2*N+1)) - dAhat_dq(:,:,1))*B;
    dBhat_dq(:,1,2) = -A\(dA_dq(:,:,2)/A*(Ahat - eye(2*N+1)) - dAhat_dq(:,:,2))*B;
    %dBhat_dq(:,1,3) = -A\(dA_dq(:,:,3)/A*(Ahat - eye(2*N+1)) - dAhat_dq(:,:,3))*B;
    dBhat_dq(:,1,3) = -A\(dA_dq(:,:,3)/A*(Ahat - eye(2*N+1)) - dAhat_dq(:,:,3))*B...
        +(Ahat-eye(2*N+1))*(A\[zeros(N+1,1);-2*N/ell^2;zeros(N-1,1)]);
end