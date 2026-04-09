function [out1, out2] = kfUpdate(arg1, arg2, arg3, arg4, arg5)
% kfUpdate  One-step linear Kalman Filter measurement update.
%
% Shorthand form — Measurement struct carries the model (C and H inferred):
%   est_state          = kfUpdate(s_pred, msmt)
%   [x_est, P_est]     = kfUpdate(x_pred, P_pred, msmt)
%
% State-struct form — update a predicted State struct:
%   est_state = kfUpdate(s_pred, zeta_or_msmt, C, H)
%
% Explicit form — operate directly on a state vector and covariance:
%   [x_est, P_est] = kfUpdate(x_pred, P_pred, zeta_or_msmt, C, H)
%
% INPUTS (shorthand forms)
%   msmt          Measurement struct from makeMeasurement whose .msmt_model
%                 field is non-empty.  C is taken from msmt_model.R and H is
%                 evaluated by calling msmt_model.h_fun (struct form) or
%                 msmt_model.h_fun_raw (explicit form) at the predicted state.
%
% INPUTS (state-struct form)
%   s_pred        Predicted State struct from predictState / kfPredict
%   zeta_or_msmt  Measurement vector (num_msmt x 1), or a Measurement struct
%                 from makeMeasurement (zeta is extracted from .zeta)
%   C             Measurement noise covariance (num_msmt x num_msmt)
%   H             Linear measurement matrix (num_msmt x n)
%
% INPUTS (explicit form)
%   x_pred        Predicted state vector (n x 1)
%   P_pred        Predicted state error covariance (n x n)
%   zeta_or_msmt  Measurement vector (num_msmt x 1), or a Measurement struct
%   C             Measurement noise covariance (num_msmt x num_msmt)
%   H             Linear measurement matrix (num_msmt x n)
%
% OUTPUTS (state-struct / shorthand-struct forms)
%   est_state   Updated State struct (same state_space and time as s_pred)
%
% OUTPUTS (explicit / shorthand-explicit forms)
%   x_est       Updated state vector (n x 1)
%   P_est       Updated state error covariance (n x n)
%
% Nicholas O'Donoughue
% June 2025

if isstruct(arg1)
    % State-struct branch: first arg is a State struct
    s_pred     = arg1;
    x_pred     = s_pred.state;
    P_pred     = s_pred.covar;
    use_struct = true;

    if nargin == 2 && isstruct(arg2) && isfield(arg2, 'msmt_model')
        % Shorthand: (s_pred, msmt) — infer C and H from msmt.msmt_model
        zeta = arg2.zeta;
        C    = arg2.msmt_model.R;
        H    = arg2.msmt_model.h_fun(s_pred);
    else
        % Full state-struct form: (s_pred, zeta_or_msmt, C, H)
        zeta_or_msmt = arg2;
        C = arg3;
        H = arg4;
        if isstruct(zeta_or_msmt)
            zeta = zeta_or_msmt.zeta;
        else
            zeta = zeta_or_msmt(:);
        end
    end
else
    % Explicit branch: first arg is a raw state vector
    x_pred     = arg1;
    P_pred     = arg2;
    use_struct = false;

    if nargin == 3 && isstruct(arg3) && isfield(arg3, 'msmt_model')
        % Shorthand: (x_pred, P_pred, msmt) — infer C and H from msmt.msmt_model
        zeta = arg3.zeta;
        C    = arg3.msmt_model.R;
        H    = arg3.msmt_model.h_fun_raw(x_pred);
    else
        % Full explicit form: (x_pred, P_pred, zeta_or_msmt, C, H)
        zeta_or_msmt = arg3;
        C = arg4;
        H = arg5;
        if isstruct(zeta_or_msmt)
            zeta = zeta_or_msmt.zeta;
        else
            zeta = zeta_or_msmt(:);
        end
    end
end

% h_fun may return a 3-D array; squeeze and transpose to (num_msmt x num_states).
if ndims(H) == 3
    H = H(:, :, 1)';
end

% KF update
z_hat = H * x_pred;
y     = zeta(:) - z_hat(:);
S     = H * P_pred * H' + C;
K     = P_pred * H' / S;
x_est = x_pred + K * y;
P_est = (eye(size(P_pred)) - K * H) * P_pred;

if use_struct
    out1 = tracker.makeState(s_pred.state_space, s_pred.time, x_est, P_est);
    out2 = [];
else
    out1 = x_est;
    out2 = P_est;
end
