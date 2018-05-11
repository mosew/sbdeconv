function KN = build_KN(a0,a1)
    global scrKN scrKNm N
    KN = blkdiag(diag([a0,zeros(1,N)])+scrKN, scrKNm);
    KN(1,N+1) = KN(1,N+1)+a1;
    KN(1,end) = KN(1,end)+0.1;
end