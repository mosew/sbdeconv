ff = figure;
for i = 1:4
    subplot(2,2,i)
    plot(t,ys(i,:),'r.')
    hold on
    for k = 1:3
        N = 2^(k+2);
        common
        uu = b{k,i}.full_deconvolved_BrAC;
        tt = t(1:(32/N):end);
        plot(tt,delaymodel(a0_stars(k,i),ell_stars(k,i),uu,tt));
    end
    legend('actual','N=8','N=16','N=32')
end