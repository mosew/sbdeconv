function AN = build_A_qq(q1M)
    global M_state L N
    Kq = build_Kq(q1M);
    AN=-M_state\(L+Kq);
    %AN = -M_state\(L+R+Kq); % different boundary condition here leads to different AN.
    assert(all(size(AN)==[N+1,N+1]));

end