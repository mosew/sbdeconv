function MN = build_MN(ell)
    global N scrMNdell scrMNmdell
    MN = blkdiag(diag([1,zeros(1,N)])+scrMNdell*ell, scrMNmdell*ell);
end