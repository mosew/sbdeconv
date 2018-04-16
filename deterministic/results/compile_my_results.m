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
global tau u_total n PAD P NUM_EPISODES

globals

test_us = u_total(:,1:end-PAD);

%% Define cell array to hold results.
fprintf('Creating empty cell array\n')
b = cell(3,NUM_EPISODES);
fprintf('Done creating empty cell array\n')
numRuns=numel(b);
thisRun=0;
rtTot=0;

%% Run tests
for i = 1:NUM_EPISODES
    
    test=i;

    % For testing paradigms 1 through 5
    % paradigm 1: testing on i, training on i+2 (wraparound)
    % paradigm 2: testing on i, training on i+1 : i+4 (wraparound)
    % paradigm 3: training on all except test episode
    
    for para = 1:3
        if para == 1
            training = rem(i+1,NUM_EPISODES)+1;
        end
        if para==2
            if i<=(NUM_EPISODES-5)
                training = (i+1):i+4;
            else
                training = [(i+1):min(i+4,NUM_EPISODES),1:(i+4-NUM_EPISODES)];
            end
        end
        if para == 3
            training = [1:(i-1),(i+1):NUM_EPISODES];
        end
        
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
        %u_star = u_star(1:(end-PAD));
        test_u = test_us(test);
        [peak_est, peaktime_est] = max(u_star);
        [peak_act, peaktime_act] = max(test_u);

        % Define and collect struct of collected data
        fprintf('Filling in cell %i,%i\n\n',para,test);
        b{para,test} = struct('tau',{tau},'P',{P},...
                                 'training_episodes',{training},'test_episode',{test},...
                                 'trained_parameters',{[q1_star,q2_star]},...
                                 'trained_regularization_scale',{[lambda1_star,lambda2_star]},...
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