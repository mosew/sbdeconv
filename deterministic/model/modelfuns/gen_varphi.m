function varphi = gen_varphi(ANhat,BNhat,u)
    global N n
    varphi=zeros(N+1,n);
    % goes from t = 0 to t = tau*(n-1), where t = (j-1)*tau --> j = floor(t/tau) + 1
    for j=2:n
        varphi(:,j)=ANhat * varphi(:,j-1) + BNhat * u(j-1);
    end
end