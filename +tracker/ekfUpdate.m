function [out1, out2] = ekfUpdate(arg1, arg2, arg3, arg4, arg5, arg6)
% ekfUpdate  One-step Extended Kalman Filter measurement update.
%
% Shorthand form — Measurement struct carries the model (C, z_fun, h_fun inferred):
%   est_state = ekfUpdate(s_pred, msmt)
%
% State-struct form — update a predicted State struct:
%   est_state = ekfUpdate(s_pred, zeta_or_msmt, C, z_fun, h_fun)
%
% Explicit form — operate directly on a state vector and covariance:
%   [x_est, P_est] = ekfUpdate(x_pred, P_pred, zeta_or_msmt, C, z_fun, h_fun)
%
% INPUTS (shorthand form)
%   msmt          Measurement struct from makeMeasurement whose .msmt_model
%                 field is non-empty.  C, z_fun, and h_fun are extracted
%                 automatically from msmt.msmt_model.
%
% INPUTS (state-struct form)
%   s_pred        Predicted State struct from predictState / ekfPredict
%   zeta_or_msmt  Measurement vector (num_msmt x 1), or a Measurement struct
%                 from makeMeasurement (zeta is extracted from .zeta)
%   C             Measurement noise covariance (num_msmt x num_msmt)
%   z_fun         Measurement function handle: z_hat = z_fun(s_pred)
%                 where s_pred is the full State struct
%   h_fun         Jacobian function handle: H = h_fun(s_pred), shape (num_msmt x n)
%
% INPUTS (explicit form)
%   x_pred        Predicted state vector (n x 1)
%   P_pred        Predicted state error covariance (n x n)
%   zeta_or_msmt  Measurement vector (num_msmt x 1), or a Measurement struct
%   C             Measurement noise covariance (num_msmt x num_msmt)
%   z_fun         Measurement function handle: z_hat = z_fun(x_pred)
%                 where x_pred is the raw state vector (n x 1)
%   h_fun         Jacobian function handle: H = h_fun(x_pred), shape (num_msmt x n)
%
% OUTPUTS (state-struct / shorthand-struct forms)
%   est_state   Updated State struct (same state_space and time as s_pred)
%
% OUTPUTS (explicit / shorthand-explicit forms)
%   x_est       Updated state vector (n x 1)
%   P_est       Updated state error covariance (n x n)
%
% Note: in the state-struct form z_fun and h_fun receive the full State struct;
% in the explicit form they receive the raw state vector (n x 1).
%
% Nicholas O'Donoughue
% June 2025

if isstruct(arg1)
    % State-struct branch: first arg is a State struct
    s_pred     = arg1;
    P_pred     = s_pred.covar;
    use_struct = true;

    if nargin == 2 && isstruct(arg2) && isfield(arg2, 'msmt_model')
        % Shorthand: (s_pred, msmt) — infer C, z_fun, h_fun from msmt.msmt_model
        zeta  = arg2.zeta;
        C     = arg2.msmt_model.R;
        z_fun = arg2.msmt_model.z_fun;
        h_fun = arg2.msmt_model.h_fun;
    else
        % Full state-struct form: (s_pred, zeta_or_msmt, C, z_fun, h_fun)
        zeta_or_msmt = arg2;
        C     = arg3;
        z_fun = arg4;
        h_fun = arg5;
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

        % Full explicit form: (x_pred, P_pred, zeta_or_msmt, C, z_fun, h_fun)
        zeta_or_msmt = arg3;
        C     = arg4;
        z_fun = arg5;
        h_fun = arg6;
        if isstruct(zeta_or_msmt)
            zeta = zeta_or_msmt.zeta;
        else
            zeta = zeta_or_msmt(:);
        end
end

% EKF update — linearise at the predicted state
if use_struct
    z_hat  = z_fun(s_pred);
    H      = h_fun(s_pred);
    x_pred = s_pred.state;
else
    z_hat = z_fun(x_pred);
    H     = h_fun(x_pred);
end

% h_fun may return a 3-D array (num_states x num_msmt x 1) for single-source
% measurement models; squeeze and transpose to (num_msmt x num_states).
if ndims(H) == 3
    H = H(:, :, 1)';
end

y     = zeta(:) - z_hat(:);
S     = H * P_pred * H' + C;
K     = P_pred * H' / S;
x_est = x_pred + K * y;
P_est = (eye(size(P_pred)) - K * H) * P_pred;

% Enforce symmetry and positive semi-definiteness after the EKF update.
% The non-Joseph form (I-KH)*P can accumulate asymmetry and small negative
% eigenvalues, especially after large-innovation (false-alarm) updates.
% This matches Python CovarianceMatrix.ensure_positive_definite(tolerance=1e-10).
P_est = (P_est + P_est') / 2;
[V, D] = eig(P_est);
P_est  = V * diag(max(diag(D), 1e-10)) * V';
P_est  = (P_est + P_est') / 2;

if use_struct
    out1 = tracker.makeState(s_pred.state_space, s_pred.time, x_est, P_est);
    out2 = [];
else
    out1 = x_est;
    out2 = P_est;
end
