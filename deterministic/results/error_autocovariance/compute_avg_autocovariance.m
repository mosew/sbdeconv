%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executable.
% 
% Used after optimization to compute mean autocovariance of process noise,
% which will be used in stochastic deconvolution algorithm.
%

load('../../data/1splhr_results.mat','actual_errors')
para = 3;

dy = zeros(9,413);
figure;hold on
for i = 1:9
    subplot(3,3,i)
    d = actual_errors{3,i};
    plot(d)
    dy(i,:) = xcov(d);
end
ACVy = mean(dy(:,length(d):end),1);
ACVy(97:end)=0;
csvwrite('acvy.csv',ACVy);
figure
plot(ACVy,'.')