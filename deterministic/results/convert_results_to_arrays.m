% Converts test results struct into MATLAB arrays.
clear
load('data\temp_results.mat','b')
load('data\preprocessed.mat','t_out')

bsize = size(b);
L2_errors = zeros(bsize);
Linf_errors = zeros(bsize);
AUC_abs_errors = zeros(bsize);
peak_time_abs_errors = zeros(bsize);
peak_height_abs_errors = zeros(bsize);
initial_q = cell(bsize);
lambdas = cell(bsize);
actual_errors = cell(bsize);
full_deconvolved_BrACs = cell(bsize);
trained_parameters = cell(bsize);
q1s = cell(bsize);
q2s = zeros(bsize);


for para = 1:bsize(1)
    for test = 1:bsize(2)                    
        s=b{para,test};
        L2_errors(para,test) = s.L2_error;
        Linf_errors(para,test) = s.Linf_error;
        AUC_abs_errors(para,test) = s.AUC_abs_error;
        initial_q{para,test}=s.initial_q;
        lambdas{para,test}=s.lambdas;
        peak_time_abs_errors(para,test) = s.peak_time_abs_error;
        peak_height_abs_errors(para,test) = s.peak_height_abs_error;
        actual_errors{para,test} = s.actual_error;
        full_deconvolved_BrACs{para,test} = s.full_deconvolved_BrAC;
        trained_parameters{para,test} = s.trained_parameters; 
        q1s{para,test} = s.trained_parameters(1:end-1);
        q2s(para,test) = s.trained_parameters(end);
        %subplot(kkkk,llll,test)
        %hold on
        %plot(s.full_deconvolved_BrAC)
        %if i == 1
        %    legend('u','y','1 training','3 training','8 training')
        %end
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
peak_time_abs_meanerr = mean(peak_time_abs_errors,2);
peak_height_abs_meanerr = mean(peak_height_abs_errors,2);

L2_MSE_sd = std(L2_errors,0,2);
logL2_sd = std(logL2_errors,0,2);
Linf_sd = std(Linf_errors,0,2);
AUC_MSE_sd = std(AUC_abs_errors.^2,0,2);
peak_time_MSE_sd = std(peak_time_abs_errors.^2,0,2);
peak_height_MSE_sd = std(peak_height_abs_errors.^2,0,2);

global resultspath
save(resultspath)
%delete data\temp_results.mat

load('data/preprocessed.mat','u_total','y_total','kkkk','llll','q')
h = figure;
for i = 1:bsize(2)
    subplot(kkkk,llll,i)
    plot(t_out,u_total(i,:),'k+')
    hold on
    plot(t_out,y_total(i,:),'r.')
    for k = 1:bsize(1)
        uu = b{k,i}.full_deconvolved_BrAC;
        plot(t_out,uu)
%         hold on
%         plot(t_out,mymodel([q1s{k,i},q2s(k,i)],uu),'gx');
        if k == bsize(1)
            xlabel('hours')
            ylabel('$10 \times$ alcohol concentration (\%)')
        end
        %ttext = ['q_1=',num2str(q(i,1)),'x+',num2str(q(i,2)),', q_2=',num2str(q(i,3))];
        %ttext = ['q_1=',num2str(q(i,1)),', q_2=',num2str(q(i,2))];
        %title(ttext,'Interpreter','tex');
        ylim([0,1.3*max(u_total(i,:))]);
    end
end
h.PaperPositionMode = 'auto';
set(0,'defaultTextInterpreter','latex')
suptitle('Estimated BAC from TAC, artificially generated TAC, N=8, $q_1=1$, $q_2=.3$, $\lambda_1=.0002$, $\lambda_2=.00001$')
%legend('BrAC','Simulated TAC','M=1', 'M=2', 'M=4', 'M=8')
%legend('BrAC','Simulated TAC','2 spl/hr', '3 spl/hr', '4 spl/hr', '6 spl/hr')
%legend('BrAC','Simulated TAC','N=4','N=8','N=16','N=32','N=64')
legend('BrAC','Simulated TAC','1 training','3 training','8 training')
