function Kq = build_Kq(q1)
    % This matrix gives [int_0^1 q^M_1(x) psi'^N_i(x) psi'N_j(x) ]_{ij}
    global N
    Kq = q1*N*(-diag(ones(1,N),-1) + [1,2*diag(ones(1,N-1)),1] - diag(ones(1,N),1));
end