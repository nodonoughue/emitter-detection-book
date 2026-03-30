function [x_pred, P_pred] = kfPredict(x_est, P_est, Q, F, u)
% kfPredict  One-step linear Kalman Filter prediction.
%
% [x_pred, P_pred] = kfPredict(x_est, P_est, Q, F)
% [x_pred, P_pred] = kfPredict(x_est, P_est, Q, F, u)
%
% INPUTS
%   x_est   State estimate vector (n x 1)
%   P_est   State error covariance (n x n)
%   Q       Process noise covariance (n x n)
%   F       Transition matrix (n x n)
%   u       Control input vector (n x 1), optional.  When provided, added
%           to the predicted mean: x_pred = F*x + u.  Does not affect the
%           covariance propagation.  Pass [] or omit to skip.
%
% OUTPUTS
%   x_pred  Predicted state vector (n x 1)
%   P_pred  Predicted state error covariance (n x n)

if nargin < 5 || isempty(u)
    u = zeros(size(x_est));
end

x_pred = F*x_est + u;

P_pred = F * P_est * F' + Q;