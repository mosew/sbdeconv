function BN = build_B(q2)
    global M_state N
    BN = M_state\[q2;zeros(N,1)];
    assert(all(size(BN)==[N+1,1]));
end