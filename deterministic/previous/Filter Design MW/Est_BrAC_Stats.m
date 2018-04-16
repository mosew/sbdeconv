function [Max_BrAC,Time_Max_BrAC,Area_BrAC_Curve] = Est_BrAC_Stats(Est_BrAC_TAC)

% This routine assume that the input is exactly as it comes out of
% BrAC_Est_0_G_1_FD, in other words that the times have NOT been
% shifted to account for the TAC latency. 
%   

[Max_BrAC,I] = max(Est_BrAC_TAC(:,2));
Time_Max_BrAC = Est_BrAC_TAC(I,1);
Area_BrAC_Curve = trapz(Est_BrAC_TAC(:,1),Est_BrAC_TAC(:,2));
end

