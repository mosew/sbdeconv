function dBNhat_dq = build_dBNhat_dq(AN,ANhat,dANhat_dq,BN)
    global N dAN_dq dBN_dq
    dBNhat_dq = zeros(N+1,1,2);
    
    dBNhat_dq(:,1,1) = -AN\(dAN_dq(:,:,1)/AN*(ANhat - eye(N+1)) - dANhat_dq(:,:,1))*BN;
    dBNhat_dq(:,1,2) =AN\(ANhat-eye(N+1))*dBN_dq(:,:,2);
end