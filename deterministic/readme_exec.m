% % readme_exec.m


%% Semi-blind deconvolution in linear time-invariant systems.


% x(j+1) = Ax(j) + Bu(j)
% y(j)   = Cx(j)            (*)
% x(0)   = 0
% 
% for j = 1,...,n
% 
% where each x(j) is an N+1-vector representing a state,
% u(j) is the average value of the signal u on [(j-1)*tau, j*tau],
% and the operators A,B,C depend on a parameter q = (q1,q2).
 

% Write
%
% u = phi*c
%
% where c=[c(1) .... c(P)]' is a P-vector and
% phi is a rank P, n x P matrix.
%       phi represents the AVE approximations of the linear splines on a
%       P-mesh of the interval [0,(n-1)*tau].
% dphi is the AVE approximations of the derivatives of the linear splines.

% Then
%
% y = L(q,u) = L(q,phi*c) := [  sum_{k=0}^{j-1}  CA^{j-k-1}B*(phi*c)  ]_{j=1}^n



% MINIMIZATION PROBLEM:

% We observe paired training data
% (u(i,j),y(i,j)) at times j=1,...,n
% for separate episodes i=1,..,m.

% Given this paired training data and a potentially unseen
% output signal y=[y(1),...,y(n)], determine the inputs u(j).

% J(q(1),q(2),c(1),...c(p),lambda1,lambda2) = 
%     1/m * sum_{i=1}^m sum_{j=1}^n (L(q,u(i,j)) - y(i,j))^2 +
%                       sum_{j=1}^n (L(q, phi*c)  - y0(j))^2
%                + lambda1* ||phi*c||^2 
%                + lambda2* ||dphi*c||^2.


% Entire folder and subfolders should be in MATLAB path.
addpath(genpath("."));

%% data/ 
% 	contains raw and processed data

%%  preprocessing/ 
% 	MATLAB scripts for preprocessing BrAC and TAC data .mat
global DATA_FILEPATH
DATA_FILEPATH = 'data\ZD_Data_5122.mat';
%nSPLHR = 6; % Approx number of splines to use per hour, assuming time is measured in hours
%preprocess
%generate_artificial_data

%% model/
% 	Contains definitions of system operators and cost functional for optimization problem.
globals

%% optimization/
% 	Contains optimization function

%% results/
% 	Calls function in optimization/ using model/ on data/preprocessed.mat and compiles error measures

global resultspath
resultspath = 'data/3splhr_results_test2_person.mat';
compile_my_results
convert_results_to_arrays
v = [L2_MSE,L2_MSE_sd,Linf_error_means,Linf_sd,AUC_MSE,AUC_MSE_sd,peak_time_MSE,peak_time_MSE_sd,peak_height_MSE,peak_height_MSE_sd]
q1v = [q(:,1)';q1s]
q2v = [q(:,2)';q2s]
save('data/3splhr_results_test2_person.mat','v','q1v','q2v');

%% previous/
% 	Contains methods of Rosen et al. for this problem, which optimizes over q, then u.
%load('rosenerrors.mat','rosenerrors')
