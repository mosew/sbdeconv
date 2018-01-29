% Use past algorithm to deconvolve BrAC from TAC for comparison

bt = cell(1,9);
errors = cell(1,9);
load('../../data/preprocessed.mat','u_total','y_total','n')

for i = 1:9
    [~,bt{i},peak_est,peaktime_est,~] = rosencode(i,rem(i+1,9)+1);

    ii = bt{i}(:,1)>=0;
    bt{i} = bt{i}(ii,:);
    
    test_u = fit_BrACTAC_to_splines({bt{i}(:,1)},{bt{i}(:,1)},{bt{i}(:,3)},{bt{i}(:,2)},1/60);
    
    
    bt{i}(:,2:3) = bt{i}(:,2:3)/100;
    
    u_star = bt{i}(1:5:length(test_u),2)';
    test_u = test_u(1:5:end);
    [peak_act, peaktime_act] = max(test_u);
    
    
        
    errors{i} = struct('tau',{tau},...
                       'training_episodes',{i+1},'test_episode',{i},...
                       'full_deconvolved_BrAC',{bt{i}(:,2)},...
                       'actual_error',{u_star-test_u},...
                       'L2_error',{sum((u_star-test_u).^2)},...
                       'Linf_error',{max(abs(u_star-test_u))},...
                       'AUC_sq_error',{(sum(u_star)-sum(test_u)).^2},...
                       'peak_time_sq_error',{(tau*(peaktime_est-peaktime_act)).^2},...
                       'peak_height_sq_error',{(peak_est-peak_act).^2});
                   %,...'trained_parameters',{[q1_star,q2_star]}
    subplot(3,3,i)
    hold on
    plot(test_u)
    plot(bt{i}(:,2))
end
