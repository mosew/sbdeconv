
function [BrAC_test,Est_BrAC_TAC,Peak_BrAC,Time_Peak_BrAC,Area_BRAC] = rosencode(training,test)
    load('AL_Data_5122.mat');
    Drink_Diary_Data_5122_A = Drink_Diary_Data_5122_A([1,2,4:8,10,11]);

    %Train filter on a single episode only
    BrAC_1 = [Drink_Diary_Data_5122_A(training).t_BrAC',Drink_Diary_Data_5122_A(training).data_BrAC];
    TAC_1 = [Drink_Diary_Data_5122_A(training).t_TAC',Drink_Diary_Data_5122_A(training).data_TAC];
    [r1_r2_h_1]=BrAC_Estimator_Filter_Design(BrAC_1,TAC_1);

    BrAC_test = [Drink_Diary_Data_5122_A(test).t_BrAC',Drink_Diary_Data_5122_A(test).data_BrAC];
    TAC_test = [Drink_Diary_Data_5122_A(test).t_TAC',Drink_Diary_Data_5122_A(test).data_TAC];
    [Est_BrAC_TAC,Peak_BrAC,Time_Peak_BrAC,Area_BRAC] = BrAC_Est_0_G_1_FD(TAC_test,r1_r2_h_1);
end

