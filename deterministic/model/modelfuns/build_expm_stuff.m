function [ANhat,dANhat_dq] = build_expm_stuff(AN)
    global dA_dq tau N
    ANhat_and_dANhatdqM = zeros(2*(N+1),2*(N+1),2);
    mat = [tau*AN,zeros(N+1);zeros(N+1),tau*AN];
    for k=1:2
        mat(1:(N+1),(N+2):end) = dA_dq(:,:,k);
        ANhat_and_dANhatdqM(:,:,k) = expm(mat);
    end
    m = ANhat_and_dANhatdqM(:,:,1)*[zeros(N+1);eye(N+1)];
    ANhat = m((N+2):end,:);
    dANhat_dq = zeros(N+1,N+1,2);
    dANhat_dq(:,:,1) = m(1:(N+1),:);
end