function mymodel = mymodel(q,u)
    global NUM_EPISODES N
    m = NUM_EPISODES;
    %m=1;
    n = size(u,2);
    state = zeros(N+1,n,m);
    mymodel = zeros(m,n);

    for j = 1:m
        q1 = q(j,1);
        q2 = q(j,2);
        A = build_A(q1);
        [Ahat,~]=build_expm_stuff(A);
        Bhat = build_Bhat(A,Ahat,build_B(q2));
        state(:,1,j) = Ahat * state(:,1,j) + Bhat * u(j,1);
        for i = 2:n
            state(:,i,j) = Ahat * state(:,i-1,j) + Bhat*u(j,i-1);
            mymodel(j,i) = state(N+1,i,j);
        end
    end
end