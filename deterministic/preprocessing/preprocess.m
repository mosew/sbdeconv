%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess.m
% 
% Executable file which preprocesses .mat file of input/output data.

% The input/output data is assumed to occupy 4 cell arrays:
% 't_TAC_5122','t_BrAC_5122','data_TAC_5122','data_BrAC_5122'
% Time is assumed measured in HOURS elapsed since episode beginning.

% %% Set preprocessing constants

global DATA_FILEPATH PAD nSPLHR
% If executing this file outside readme_exec.m, uncomment these.
% matlabpath(pathdef);
% DATA_FILEPATH = 'data\ZD_Data_5122.mat';

nSPLHR = 1; % Approx number of splines to use per hour
TAU_OUT = 5/60; % Number of hours between timesteps in spline fit
PAD = 11; % Number of zeros to add to the beginning and end of each sample (equivalent to 5*PAD minutes)

% Later the PAD entries at the end of each episode will be dropped to mitigate overshoot artifacts at end of each episode.

%% Load original data
load(DATA_FILEPATH,...
     't_TAC_5122','t_BrAC_5122','data_TAC_5122','data_BrAC_5122');
 
% %% Fix an entry from spreadsheet
% % From the spreadsheet, this was read as a 32 but it was thought unlikely.
% t_BrAC_5122{8} = [0,0.5,t_BrAC_5122{8}(2:end)];
% data_BrAC_5122{8} = [0,32,data_BrAC_5122{8}(2:end)']';

%% Exclude episodes 3 and 9
% Episode 3 the TAC signal never returns to 0
% Episode 9 Susan takes off her sensor during the drinking episode    
dontuse = [3,9];
toUseIndex = true(1,11);
toUseIndex(dontuse)=false;
t_TAC_5122 = t_TAC_5122(toUseIndex);
t_BrAC_5122 = t_BrAC_5122(toUseIndex);
data_TAC_5122 = data_TAC_5122(toUseIndex);
data_BrAC_5122 = data_BrAC_5122(toUseIndex);

global NUM_EPISODES
NUM_EPISODES = 9;

%% Pads the beginning and end of each training episode (sample) with PAD/12 hours of zeros
% AND normalizes data
for i = 1:length(data_BrAC_5122)
    t_TAC_5122{i} = [(-PAD:-1)/12, t_TAC_5122{i}, t_TAC_5122{i}(end)+(1:PAD)/12] + PAD/12;
    t_BrAC_5122{i} = [(-PAD:-1)/12,t_BrAC_5122{i}, t_BrAC_5122{i}(end)+(1:PAD)/12] + PAD/12;
    data_BrAC_5122{i} = [zeros(1,PAD),data_BrAC_5122{i}', zeros(1,PAD)]/100;
    data_TAC_5122{i} = [zeros(1,PAD),data_TAC_5122{i}',zeros(1,PAD)]/100;
end

%% Interpolate data
% BrAC (input) data is interpolated with splines
% TAC (output) data is interpolated with linear splines, but in such a way
% as to preserve the area under the curve, because the TAC samples
% represent accumulations of molecules over a 5-minute window.

[t_out,u_total,y_total] = fit_BrACTAC_to_splines(t_TAC_5122,t_BrAC_5122,data_TAC_5122,data_BrAC_5122,TAU_OUT);
tau = TAU_OUT;
n = length(t_out);
kkkk = fix(sqrt(NUM_EPISODES));
llll = ceil(NUM_EPISODES/kkkk);

save('data\preprocessed.mat','u_total','y_total','tau','n','PAD','NUM_EPISODES','nSPLHR','kkkk','llll');

%% Check data
for i = 1:NUM_EPISODES
    subplot(kkkk,llll,i)
    plot(u_total(i,:),'--')
    hold on
    plot(y_total(i,:),'.')
end
suptitle('Preprocessed data')
legend('BrAC','TAC')

clear