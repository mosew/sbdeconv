% % readme_exec.md


%% Semi-blind deconvolution in linear time-invariant systems.

% x(j+1) = Ax(j) + Bu(j)
% y(j)   = Cx(j)            (*)
% x(0)   = 0
% 
% for j = 1,...,n
% 
% where each x(j) is an N+1-vector representing a state.
% and the operators A,B,C depend on a parameter q = (q1,q2).
% 

% Write
%
% u = phi*c
%
% where c=[c(1) .... c(P)]' is a P-vector and
% phi is a rank P, n x P matrix.
%       phi represents the AVE approximations of the linear splines on a
%       P-mesh of the interval [0,n*tau].
% dphi is the AVE approximations of the derivatives of the linear splines.

% Then
% y = L(q,u) = L(q,phi*c) := [sum_{k=0}^{j-1} CA^{j-k-1}B*(phi*c)]_{j=1}^n



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
addpath(genpath("C:\Users\mose\Dropbox\research\deterministic"));


%% data/ 
% 	contains raw and processed data

%%  preprocessing/ 
% 	MATLAB scripts for preprocessing BrAC and TAC data from Excel
% 	spreadsheet and .mat
global DATA_FILEPATH TAU_IN TAU_OUT PAD
DATA_FILEPATH = '\data\ZD_Data_5122.mat';
TAU_IN = 5/60; % Number of hours between timesteps in TAC sample
TAU_OUT = 5/60; % Number of hours between timesteps in spline fit
PAD = 11; % Number of zeros to add to the beginning and end of each sample (equivalent to 5*PAD minutes)
% Later the PAD entries at the end of each episode will be dropped to
% mitigate overshoot artifacts at end of each episode.
preprocess

%% model/
% 	Contains definitions of system operators and cost functional for optimization problem.
JN_and_dJN_globals

%% optimization/
% 	Contains optimization function

%% results/
% 	Calls function in optimization/ using model/ on data/preprocessed.mat and compiles error measures

global resultspath
resultspath = 'data/results_arrays.mat';
compile_my_results
convert_results_to_arrays

%% previous/
% 	Contains methods of Rosen et al. for this problem, which optimizes over q, then u.
 