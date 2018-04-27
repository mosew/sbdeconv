function dA_dq = build_dAN_dqM()
    global M_state N M
    dA_dq = zeros(N+1,N+1,M+2);
    % Last page is q2, all zeros.
    for k=1:M+1
        z = zeros(1,M+1);
        z(k)=1;
        dA_dq(:,:,k) = -M_state\build_Kq(z);
    end
end