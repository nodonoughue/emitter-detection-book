function [out1, out2] = kfPredict(arg1, arg2, Q, F, u)
% kfPredict  One-step linear Kalman Filter prediction.
%
% State-struct form — predict a State struct forward to new_time:
%   s_pred = kfPredict(state, new_time, Q, F)
%   s_pred = kfPredict(state, new_time, Q, F, u)
%
% Explicit form — operate directly on a state vector and covariance:
%   [x_pred, P_pred] = kfPredict(x_est, P_est, Q, F)
%   [x_pred, P_pred] = kfPredict(x_est, P_est, Q, F, u)
%
% INPUTS (state-struct form)
%   state     State struct from makeState (must have .state, .covar, .state_space)
%   new_time  Target timestamp [s] for the predicted state
%   Q         Process noise covariance (n x n)
%   F         Transition matrix (n x n)
%   u         Control input vector (n x 1), optional.  Added to the predicted
%             mean: x_pred = F*x + u.  Does not affect covariance.  Default: []
%
% INPUTS (explicit form)
%   x_est     State estimate vector (n x 1)
%   P_est     State error covariance (n x n)
%   Q         Process noise covariance (n x n)
%   F         Transition matrix (n x n)
%   u         Control input vector (n x 1), optional.  Default: []
%
% OUTPUTS (state-struct form)
%   s_pred    Predicted State struct at new_time
%
% OUTPUTS (explicit form)
%   x_pred    Predicted state vector (n x 1)
%   P_pred    Predicted state error covariance (n x n)
%
% Nicholas O'Donoughue
% June 2025

if isstruct(arg1)
    % State-struct form: (state, new_time, Q, F, u)
    state    = arg1;
    new_time = arg2;
    if nargin < 5 || isempty(u)
        u = zeros(size(state.state));
    end
    x_pred = F * state.state + u;
    P_pred = F * state.covar * F' + Q;
    out1 = tracker.makeState(state.state_space, new_time, x_pred, P_pred);
    out2 = [];
else
    % Explicit form: (x_est, P_est, Q, F, u)
    x_est = arg1;
    P_est = arg2;
    if nargin < 5 || isempty(u)
        u = zeros(size(x_est));
    end
    out1 = F * x_est + u;
    out2 = F * P_est * F' + Q;
end
