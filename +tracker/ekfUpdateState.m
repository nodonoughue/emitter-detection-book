function s_upd = ekfUpdateState(s_pred, zeta, msmt_model)
% ekfUpdateState  EKF measurement update on a State struct.
%
% s_upd = ekfUpdateState(s_pred, zeta, msmt_model)
%
% Linearises the measurement model at s_pred, computes the Kalman gain, and
% returns the updated State.  For linear measurement models the result is
% identical to a standard KF update.
%
% INPUTS
%   s_pred      Predicted State struct (from predictState)
%   zeta        Measurement vector (num_msmt x 1)
%   msmt_model  Measurement model struct from makeMeasurementModel
%
% OUTPUTS
%   s_upd   Updated State struct (same state_space and time as s_pred)
%
% Nicholas O'Donoughue
% June 2025

% Expected measurement and Jacobian at the predicted state
z_hat = msmt_model.z_fun(s_pred);     % (num_msmt x 1)
H     = msmt_model.h_fun(s_pred);     % (num_msmt x num_states)

% Handle the case where h_fun returns a 3-D array (num_states x num_msmt x 1)
% as produced by the existing makeMeasurementModel helper when num_src=1.
if ndims(H) == 3
    H = H(:, :, 1)';   % squeeze and transpose to (num_msmt x num_states)
end

P = s_pred.covar;
R = msmt_model.R;

% Innovation
y = zeta(:) - z_hat(:);

% Innovation covariance
S = H * P * H' + R;

% Kalman gain
K = P * H' / S;

% Updated state mean and covariance
n = size(P, 1);
x_upd = s_pred.state + K * y;
P_upd = (eye(n) - K * H) * P;

s_upd = tracker.makeState(s_pred.state_space, s_pred.time, x_upd, P_upd);
