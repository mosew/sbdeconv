%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Executable.
%
% Compiles error measures for my optimization routine for various
% training paradigms, i.e. numbers of training episodes.
% 
% Be sure "model" and "optimization" folders are in MATLAB path before running.

%% Set up globals
rng(323)

clear
try
    load('data\preprocessed.mat')
catch
    preprocess
    load('data\preprocessed.mat')
end
global tau u_total P NUM_EPISODES N

globals_qq

test_us = u_total;

%% Define cell array to hold results.
fprintf('Creating empty cell array\n')
b = cell(5,9);
fprintf('Done creating empty cell array\n')
numRuns=numel(b);
thisRun=0;
rtTot=0;

%% Run tests
for i = 1:NUM_EPISODES
    
    test=i;

    % TEST 3
    % For testing paradigms 1 through 5, test on i, train on i+2,
    % N=2^(para+1)
    
    training = rem(i+1,NUM_EPISODES)+1;
    for para = 1:5
        N = 2^(1+para);
        save('data\preprocessed.mat','N','-append')
        globals
        
        % Run the minimization script and time it
        tic
        fprintf('Testing on P=%i, paradigm=%i, test episode=%i \n', P, para, test)
        [q1_star,q2_star,u_star,lambda1_star,lambda2_star] = optimize_input_and_params(training,test,P);
        r=toc;

        thisRun = thisRun +1;
        fprintf('Completed %i / %i \n', thisRun,numRuns)
        rtTot = rtTot+r;
        fprintf('Average time per run=%f s\n',rtTot/thisRun);

        % Collect data from the run
        test_u = test_us(test,:);
        [peak_est, peaktime_est] = max(u_star);
        [peak_act, peaktime_act] = max(test_u);

        % Define and collect struct of collected data
        fprintf('Filling in cell %i,%i\n\n',para,test);
        b{para,test} = struct('tau',{tau},'P',{P},'M',{M},...
                                 'training_episodes',{training},'test_episode',{test},...
                                 'trained_parameters',{[q1_star,q2_star]},...
                                 'full_deconvolved_BrAC',{u_star},...
                                 'actual_error',{u_star-test_u},...
                                 'L2_error',{tau^2*norm(u_star-test_u,2)},...
                                 'Linf_error',{max(abs(u_star-test_u))},...
                                 'AUC_abs_error',{tau*abs(trapz(u_star)-trapz(test_u))},...
                                 'peak_time_abs_error',{tau*abs(peaktime_est-peaktime_act)},...
                                 'peak_height_abs_error',{abs(peak_est-peak_act)});
    end
end
fprintf('Saving test results cell array\n')
save('data\temp_results.mat','b')