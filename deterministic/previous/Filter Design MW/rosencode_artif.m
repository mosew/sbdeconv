
function [BrAC_test,Est_BrAC_TAC,Peak_BrAC,Time_Peak_BrAC,Area_BRAC] = rosencode_artif(training,test)
    %Train filter on a single episode only
    load('data/preprocessed.mat','u_total','y_total','t_out','PAD')
    BrAC_1 = [t_out',u_total(training,PAD+1:end-PAD)'];
    TAC_1 = [t_out',y_total(training,PAD+1:end-PAD)'];
    [r1_r2_h_1]=BrAC_Estimator_Filter_Design(BrAC_1,TAC_1);

    BrAC_test = [t_out',u_total(test,PAD+1:end-PAD)'/100];
    TAC_test = [t_out',y_total(test,PAD+1:(end-PAD))'/100];

    [Est_BrAC_TAC,Peak_BrAC,Time_Peak_BrAC,Area_BRAC] = BrAC_Est_0_G_1_FD(TAC_test,r1_r2_h_1);
end