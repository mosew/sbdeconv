t = 0:(T/(n-1)/(32/N)):T;

f = figure;
for i = 1:bsize(2)
    subplot(2,2,i)
    plot(t,us(i,:),'k+')
    hold on
    for k = 1:bsize(1)
        uu = b{k,i}.full_deconvolved_BrAC;
        plot(t(1:(32/2^(k+1)):end),uu)
        if k == bsize(1)
            xlabel('')
            ylabel('$10 \times$ alcohol concentration (\%)')
        end
        ylim([0,1.3*max(us(i,:))]);
    end
end
f.PaperPositionMode = 'auto';
set(0,'defaultTextInterpreter','latex')
suptitle('Deconvolved input signal, hereditary system, $\lambda_1=$, $\lambda_2=$, $\lambda_3=$')
legend('actual','estimated, 1 training', 'estimated, 2 training', 'estimated, 3 training')
%legend('BrAC','Simulated TAC','2 spl/hr', '3 spl/hr', '4 spl/hr', '6 spl/hr')
%legend('BrAC','Simulated TAC','N=4','N=8','N=16','N=32','N=64')
%legend('BrAC','TAC','1 training','3 training','8 training')

ff = figure;
for j = 1:4
    subplot(2,2,j)
    plot(t,ys(j,:),'r.')
    hold on
    for jj = 1:4
        N = 2^(jj+1);
        common
        uu = b{jj,j}.full_deconvolved_BrAC;
        tt = t(1:(32/N):end);
        plot(tt,delaymodel(a0_stars(jj,j),ell_stars(jj,j),uu,tt));
    end
    legend('actual','N=4','N=8','N=16','N=32')
end