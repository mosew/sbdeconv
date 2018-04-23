% Use past algorithm to deconvolve BrAC from TAC for comparison

global tau
load('data/preprocessed.mat','tau')
bt = cell(1,9);
rosenerrors = cell(1,9);
%tau = 5/60;

for i = 1:9
    [test_u,bt{i},peak_est,peaktime_est,auc_est] = rosencode_artif(rem(i+1,9)+1,i);    
    
    u_star = bt{i}(:,2);
    
    test_u = interp1(test_u(:,1)',test_u(:,2)',bt{i}(1,1):(5/60):bt{i}(end,1),'linear','extrap')';
    test_u = max(test_u,0);
    
    u_star = u_star(1:5:end);
        
    assert(all(size(u_star)==size(test_u)));
    
    [peak_act, peaktime_act] = max(test_u);
    peaktime_act = bt{i}(1)+(peaktime_act-1)*tau;
    
    rosenerrors{i} = struct('tau',{tau},...
                       'training_episodes',{i+1},'test_episode',{i},...
                       'full_deconvolved_BrAC',{bt{i}(:,2)},...
                       'actual_error',{u_star-test_u},...
                       'L2_error',{tau^2*norm(u_star-test_u,2)},...
                       'Linf_error',{max(abs(u_star-test_u))},...
                       'AUC_abs_error',{tau*abs(trapz(u_star)-trapz(test_u))},...
                       'peak_time_abs_error',{abs(peaktime_est-peaktime_act)},...
                       'peak_height_abs_error',{abs(peak_est-peak_act)});
                   %,...'trained_parameters',{[q1_star,q2_star]}
    subplot(3,3,i)
    plot(test_u)
    hold on
    plot(u_star)
end
suptitle('Rosen et al method')

% Converts test results struct into MATLAB arrays.

bsize = size(rosenerrors);
L2_errors = zeros(bsize);
Linf_errors = zeros(bsize);
AUC_abs_errors = zeros(bsize);
peak_time_abs_errors = zeros(bsize);
peak_height_abs_errors = zeros(bsize);
actual_errors = cell(bsize);
full_deconvolved_BrACs = cell(bsize);


for para = 1:bsize(1)
    for test = 1:bsize(2)                    
        s=rosenerrors{para,test};
        L2_errors(para,test) = s.L2_error;
        Linf_errors(para,test) = s.Linf_error;
        AUC_abs_errors(para,test) = s.AUC_abs_error;
        peak_time_abs_errors(para,test) = s.peak_time_abs_error;
        peak_height_abs_errors(para,test) = s.peak_height_abs_error;
        actual_errors{para,test} = s.actual_error;
        full_deconvolved_BrACs{para,test} = s.full_deconvolved_BrAC;
    end
end

% Calculate the mean and SD of error measures across test episodes.

L2_MSE = mean(L2_errors,2);
Linfmean = mean(Linf_errors,2);
AUCmean = mean(AUC_abs_errors,2);
ptmean = mean(peak_time_abs_errors,2);
phmean = mean(peak_height_abs_errors,2);

L2_MSE_sd = std(L2_errors,0,2);
Linf_sd = std(Linf_errors,1,2);
AUC_MSE_sd = std(AUC_abs_errors,0,2);
peak_time_MSE_sd = std(peak_time_abs_errors,0,2);
peak_height_MSE_sd = std(peak_height_abs_errors,0,2);

save('rosenerrors.mat')

