% Converts test results struct into MATLAB arrays.
clear
load('data\temp_results.mat','b')

bsize = size(b);
L2_errors = zeros(bsize);
Linf_errors = zeros(bsize);
AUC_abs_errors = zeros(bsize);
peak_time_abs_errors = zeros(bsize);
peak_height_abs_errors = zeros(bsize);
actual_errors = cell(bsize);
full_deconvolved_BrACs = cell(bsize);
trained_parameters = cell(bsize);
lambda1s = zeros(bsize);
lambda2s = zeros(bsize);


load('data/preprocessed.mat','u_total','y_total','kkkk','llll')
h = figure;
for i = 1:bsize(2)
    subplot(kkkk,llll,i)
    plot(u_total(i,:),'k--')
    hold on
    plot(y_total(i,:),'r.')
    for k = 1:bsize(1)
        plot(b{k,i}.full_deconvolved_BrAC)
    end
end

for para = 1:bsize(1)
    for test = 1:bsize(2)                    
        s=b{para,test};
        L2_errors(para,test) = s.L2_error;
        Linf_errors(para,test) = s.Linf_error;
        AUC_abs_errors(para,test) = s.AUC_abs_error;
        peak_time_abs_errors(para,test) = s.peak_time_abs_error;
        peak_height_abs_errors(para,test) = s.peak_height_abs_error;
        actual_errors{para,test} = s.actual_error;
        full_deconvolved_BrACs{para,test} = s.full_deconvolved_BrAC;
        trained_parameters{para,test} = s.trained_parameters;
        lambda1s(para,test) = s.trained_regularization_scale(1);
        lambda2s(para,test) = s.trained_regularization_scale(2);        
        
        subplot(kkkk,llll,test)
        hold on
        plot(s.full_deconvolved_BrAC)
        if i == 1
            legend('u','y','1 training','3 training','8 training')
        end
    end
end

logL2_errors = log(L2_errors);
logL2_MSE = mean(logL2_errors,2);
logAUC_sq_errors = log(AUC_abs_errors);
logAUC_MSE = mean(AUC_abs_errors,2);

% Calculate the mean and SD of error measures across test episodes.

L2_MSE = mean(L2_errors,2);
Linf_error_means = mean(Linf_errors,2);
AUC_MSE = mean(AUC_abs_errors,2);
peak_time_MSE = mean(peak_time_abs_errors,2);
peak_height_MSE = mean(peak_height_abs_errors,2);

L2_MSE_sd = std(L2_errors,0,2);
logL2_sd = std(logL2_errors,0,2);
Linf_sd = std(Linf_errors,0,2);
AUC_MSE_sd = std(AUC_abs_errors.^2,0,2);
peak_time_MSE_sd = std(peak_time_abs_errors.^2,0,2);
peak_height_MSE_sd = std(peak_height_abs_errors.^2,0,2);

global resultspath
save(resultspath)
%delete data\temp_results.mat
