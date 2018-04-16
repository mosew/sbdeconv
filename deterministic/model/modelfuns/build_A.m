function AN = build_A(q1)
    global M_state K_state L N
    AN=-M_state\(L+q1*K_state);
    %AN = -M_state\(L+R+Kq); % different boundary condition here leads to different AN.
    assert(all(size(AN)==[N+1,N+1]));
end