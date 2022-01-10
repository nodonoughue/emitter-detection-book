function [P,Pe] = steadystateError(F,H,Q,R,max_num_iterations,epsilon)
% [P,Pe] = steadystateError(F,H,Q,R,max_num_iterations)
%
% Iteratively compute the steady-state prediction error for a Kalman Filter
% based on the steady-state Ricatti equation for prediction error:
%       P = Q + F( P^{-1} + H'R^{-1}H )^{-1} F'
%
% where
%   F = state transition matrix
%   H = measurement matrix
%   Q = process noise (variance of state from one update to the next)
%   R = measurement noise (variance of each set of measurements y)
%
% This is done iteratively, until the difference between successive
% estimates of F is sufficiently small, or until the max_num_iterations is
% reached.
%
% INPUTS:
%   F       N x N state transition matrix
%   H       M x N measurement matrix
%   Q       N x N covariance of process noise
%   R       M x M covariance of measurement noise
%   max_num_iterations (optional) stopping condition for iteration
%   epsilon (optional) stopping condition for iteration
%
% OUTPUTS:
%   P       N x N steady-state Prediction error covariance
%   Pe      N x N steady-state Estimation error covariance
%
% Nicholas O'Donoughue
% 11 Nov 2021

if nargin < 6 || ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 1e-6 * norm(F,'fro');
end

if nargin < 5 || ~exist('max_num_iterations','var') || isempty(max_num_iterations)
    max_num_iterations = 100;
end

isOctave = exist('OCTAVE_VERSION','builtin')~=0;

%% Initialize P
P = eye(size(Q));
dist = inf;
iter = 0;

%% Pre-compute R inverse
HRH = (H' / R) * H;

%% Loop until convergence
while iter < max_num_iterations && abs(dist) > epsilon
    % Update estimate
    P_prev = P;
    
    % Pre-invert P_prev, using knowledge that it's a symmetric matrix
    P = Q + (F / (pinv(P_prev) + HRH)) * F';
    
    % Increment the convergence counters
    dist = norm(P-P_prev,'fro');
    iter = iter + 1;
end

%% Comput Estimation Error
Pe = inv(inv(P) + HRH);