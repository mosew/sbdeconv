function [time_out,u_total,y_total] = fit_BrACTAC_to_splines(t_TAC,t_BrAC,data_TAC,data_BrAC,tau_out)

% INPUT:
% 
% t_TAC     : cell array consisting of time vectors for all TAC episodes
% t_BrAC    :              ""                               BrAC   ""
% data_TAC  :              ""          data       ""        TAC    ""
% data_BrAC :              ""                               BrAC   ""
% tau_in    : length (in hours) of one unit in input time vectors
% tau_out   : desired output timestep (in hours)
%
% Timestep must be consistent across inputs.
% Each episode must start at time index 0 (NOT 1)

% OUTPUT
% time_out  : vector of time indices used to sample splines, with timestep tau_out
% u_total   : ""                data                BrAC    ""   , fit to spline and sampled uniformly with timestep tau
% y_total   : ""                                     TAC    ""


% Total number of episodes
m_total = length(t_TAC);

% Check to make sure time index of every episode starts at t=0
for i=1:m_total
    assert(t_TAC{i}(1)==0)
    assert(t_BrAC{i}(1)==0)
end

% Find maximum episode length (in tau_in timesteps)
max_n_in=0;
T=0;
for i=1:length(t_TAC)
    max_n_in = max(max(length(t_TAC{i}),length(t_BrAC{i})),max_n_in);
    T = max(T,max(max(t_TAC{i})));
end

% Number of timesteps in output vector
n_out = floor( T / tau_out);

time_out = (0:(n_out-1))*tau_out;

u_total=zeros(m_total,n_out);
y_total=zeros(m_total,n_out);

% Generate, evaluate, restrict>=0 splines.
for i = 1:m_total
    y_total(i,:)=max(interp1(t_TAC{i},data_TAC{i},time_out,'linear','extrap'),0);
    u_total(i,:)=max(interp1(t_BrAC{i},data_BrAC{i},time_out,'spline','extrap'),0);
end