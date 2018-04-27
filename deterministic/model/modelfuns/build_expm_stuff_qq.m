function [ANhat,dANhat_dqM] = build_expm_stuff_qq(AN)
    global dA_dq tau N M
    ANhat_and_dANhatdqM = zeros(2*(N+1),2*(N+1),M+2);
    mat = [tau*AN,zeros(N+1);zeros(N+1),tau*AN];
    dANhat_dqM = zeros(N+1,N+1,M+2);
    for k=1:(M+2)
        mat(1:(N+1),(N+2):end) = dA_dq(:,:,k);
        ANhat_and_dANhatdqM(:,:,k) = expm(mat);
        m = ANhat_and_dANhatdqM(:,:,k)*[zeros(N+1);eye(N+1)];
        dANhat_dqM(:,:,k) = m(1:(N+1),:);
        ANhat = m((N+2):end,:);
    end
end