function Kq = build_Kq(q1M)
    global si_diag si_offdiag N
    for k = 1:(M+1)
        Kq(:,:,k) = N*(diag(-q1M.*si_offdiag(k,:),-1)+diag(q1M.*si_diag(k,:))+diag(-q1M.*si_offdiag(k,:),1));
    end
