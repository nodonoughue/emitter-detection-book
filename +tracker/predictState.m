function s_pred = predictState(s, new_time, motion_model)
% predictState  Predict a State struct forward to new_time.
%
% s_pred = predictState(s, new_time, motion_model)
%
% INPUTS
%   s             State struct from makeState (must have .covar set)
%   new_time      Target timestamp [s]
%   motion_model  Motion model struct from makeMotionModel
%
% OUTPUTS
%   s_pred   New State struct at new_time, with predicted state vector and covariance
%
% Nicholas O'Donoughue
% June 2025

dt = new_time - s.time;

if dt == 0
    s_pred = s;
    return;
end

Q = motion_model.q_fun(dt);

if motion_model.is_linear
    F = motion_model.f_fun(dt);
    u = [];
    if isfield(motion_model, 'b_fun') && ~isempty(motion_model.b_fun)
        u = motion_model.b_fun(dt);
    end
    [x_pred, P_pred] = tracker.kfPredict(s.state, s.covar, Q, F, u);
else
    % Nonlinear EKF predict (e.g. Constant Turn)
    f = @(x) motion_model.ct_fun(x, dt);
    g = @(x) motion_model.jacobian_fun(x, dt);
    [x_pred, P_pred] = tracker.ekfPredict(s.state, s.covar, Q, f, g);
end

s_pred = tracker.makeState(s.state_space, new_time, x_pred, P_pred);
