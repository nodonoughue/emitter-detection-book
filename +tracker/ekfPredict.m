function [out1, out2] = ekfPredict(arg1, arg2, Q, f_fun, g_fun)
% ekfPredict  One-step Extended Kalman Filter prediction.
%
% State-struct form — predict a State struct forward to new_time:
%   s_pred = ekfPredict(state, new_time, Q, f_fun, g_fun)
%
% Explicit form — operate directly on a state vector and covariance:
%   [x_pred, P_pred] = ekfPredict(x_est, P_est, Q, f_fun, g_fun)
%
% INPUTS (state-struct form)
%   state     State struct from makeState (must have .state, .covar, .state_space)
%   new_time  Target timestamp [s] for the predicted state
%   Q         Process noise covariance (n x n)
%   f_fun     Nonlinear transition function handle: x_pred = f_fun(x)
%   g_fun     Jacobian function handle: F = g_fun(x), shape (n x n)
%
% INPUTS (explicit form)
%   x_est     State estimate vector (n x 1)
%   P_est     State error covariance (n x n)
%   Q         Process noise covariance (n x n)
%   f_fun     Nonlinear transition function handle: x_pred = f_fun(x)
%   g_fun     Jacobian function handle: F = g_fun(x), shape (n x n)
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
    % State-struct form: (state, new_time, Q, f_fun, g_fun)
    state    = arg1;
    new_time = arg2;
    x_pred = f_fun(state.state);
    F      = g_fun(state.state);
    P_pred = F * state.covar * F' + Q;
    out1 = tracker.makeState(state.state_space, new_time, x_pred, P_pred);
    out2 = [];
else
    % Explicit form: (x_est, P_est, Q, f_fun, g_fun)
    x_est = arg1;
    P_est = arg2;
    x_pred = f_fun(x_est);
    F      = g_fun(x_est);
    out1 = x_pred;
    out2 = F * P_est * F' + Q;
end
