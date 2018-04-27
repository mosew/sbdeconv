function Kq = build_Kq(q1M)
    global si_diag si_offdiag N M
    Kq = zeros(N+1,N+1);
    for k = 1:(M+1)
        Kq = diag(-q1M*si_offdiag,-1)+diag(q1M*si_diag)+diag(-q1M*si_offdiag,1);
    end
end