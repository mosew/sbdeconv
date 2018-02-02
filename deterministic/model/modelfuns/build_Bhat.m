function Bhat = build_Bhat(A,Ahat,B)
    global N
    Bhat = (Ahat-eye(N+1))*(A\B);
end
