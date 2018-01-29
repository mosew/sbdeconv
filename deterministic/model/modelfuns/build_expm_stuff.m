function [ANhat,dANhat_dqM] = build_expm_stuff(AN)
    global dAN_dqM tau N
    ANhat_and_dANhatdqM = zeros(2*(N+1),2*(N+1),2);
    mat = [tau*AN,zeros(N+1);zeros(N+1),tau*AN];
    for k=1:2
        mat(1:(N+1),(N+2):end) = dAN_dqM(:,:,k);
        ANhat_and_dANhatdqM(:,:,k) = expm(mat);
    end
    m = ANhat_and_dANhatdqM(:,:,1)*[zeros(N+1);eye(N+1)];
    ANhat = m((N+2):end,:);
    dANhat_dqM = m(1:(N+1),:);
end