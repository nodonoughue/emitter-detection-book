function [P_pred, P_est] = modelKalmanStepError(x_pred, P_pred, F, H, Q, R, update_type, msmt_type, f_jacob, h_jacob)
% [P_pred, P_est] = modelKalmanStepError(x_pred, P_pred, F, H, Q, R)
%
%   Simplified Kalman Filter loop iteration; assumes a strict KF
%   implementation (linear system), or an extended KF.
%
%   Accepts a predicted state and covariance error from the previous time
%   step, and a set of transition, measurement, process noise, and 
%   measurement noise matrices that describe the system.
%
%   Models the impact of a single time step update on the predicted and
%   estimated state error covariance (P_pred, P_est).
%
%   Inputs:
%       x_pred      N x 1 predicted target state from the prior time step
%       P_pred      N x N predicted target state covariance, from the
%                   prior time step
%       F           N x N state transition matrix (or function handle for
%                   the Jacobian of F(x) is update_type is 'kf')
%       H           M x N measurement matrix (or function handle for the
%                   Jacobian of H(x) if msmt_type is 'ekf')
%                   generate it)
%       Q           N x N process noise covariance matrix
%       R           M x M measurement noise covariance matrix
%       update_type (Optional) string indicating KF type {kf, ekf} for the
%                   update (state transition) model. Default=kf
%       msmt_type   (Optional) string indicating KF type {kf, ekf} for the
%                   measurement model. Default=kf
%
% Outputs:
%       P_pred      N x N predicted target state covariance matrix
%       P_est       N x N estimated target state covariance matrix
%
% Nicholas O'Donoughue
% 4 Jan 2022

%% Parse Inputs

if nargin < 7 || isempty(update_type)
    update_type = 'kf';
end
if nargin < 8 || isempty(msmt_type)
    msmt_type = 'kf';
end

N = size(x_pred,1); % Number of states

assert(size(P_pred,1)==N && size(P_pred,2)==N,'Dimension mismatch between target state and error covariance.');
assert(strcmpi(update_type,'ekf')==1 || (size(F,1)==N && size(F,2)==N),'Dimension mismatch for state transition matrix.');
assert(size(Q,1)==N && size(Q,2)==N,'Dimension mismsatch between target and process noise covariance.');

if strcmpi(update_type,'ekf')
    assert(nargin > 8 && ~isempty(f_jacob), 'If update type is EKF, then f_jacob is required.');
end
if strcmpi(msmt_type,'ekf')
    assert(nargin > 9 && ~isempty(h_jacob), 'If measurement type is EKF, then h_jacob is required.');
end

%% State Estimation Step
% Compute the Kalman Gain and do update the current state estimate
% based on the provided measurement

% Parse KF or EKF model type
switch lower(msmt_type)
    case 'kf'
        % H is a matrix
        H_lin = H;
    case 'ekf'
        % H is a function handle, use the jacobian to compute a
        % linearized H
        H_lin = H(x_pred);
    otherwise
        error('Unrecognized measurement type. Must be ''kf'' or ''ekf''.');
end

% Compute error covariance for residual
S = H_lin*P_pred*H_lin' + R;

% Optimal Kalman Gain
K = P_pred*H_lin'/S;

% Updated State Estimate Covariance
P_est = (eye(N)-K*H)*P_pred;


%% State Prediction Step

% Parse KF or EKF model type
switch lower(update_type)
    case 'kf'
        % F is a matrix
        F_lin = F;
    case 'ekf'
        % F is a function handle, use it to compute a linearized F
        F_lin = F(x_pred);
    otherwise
        error('Unrecognized update type. Must be ''kf'' or ''ekf''.');
end

% Propagate forward the error covariance, and then add process noise
% covariance
P_pred = F_lin*P_est*F_lin' + Q;