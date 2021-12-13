function [x_pred, x_est, P_pred, P_est] = runKalmanStep(x_pred, P_pred, z, F, H, Q, R, update_type, msmt_type, f_jacob, h_jacob)
% [x_pred, x_est, P_pred, P_est] = runKalmanStep(x_pred, P_pred, z, F, H, Q, R)
%
%   Simplified Kalman Filter loop iteration; assumes a strict KF
%   implementation (linear system), or an extended KF.
%
%   Accepts a predicted state and covariance error from the previous time
%   step, measurements for the current step, and a set of
%   transition, measurement, process noise, and measurement noise
%   matrices that describe the system.
%
%   Ignores any problems involving data association (selecting which track
%   goes with which measurement), track initiation or deletion.
%
%   To model a missed detection, call this function with an empty
%   measurement.
%
%   Inputs:
%       x_pred      N x 1 predicted target state from the prior time step
%       P_pred      N x N predicted target state covariance, from the
%                   prior time step
%       z           M x 1 measurement vector for the current time step
%       F           N x N state transition matrix (or function handle to
%                   generate it)
%       H           M x N measurement matrix (or function handle to
%                   generate it)
%       Q           N x N process noise covariance matrix
%       R           M x M measurement noise covariance matrix
%       update_type (Optional) string indicating KF type {kf, ekf} for the
%                   update (state transition) model. Default=kf
%       msmt_type   (Optional) string indicating KF type {kf, ekf} for the
%                   measurement model. Default=kf
%       f_jacob     Function handle to generator for jacobian of state
%                   transition matrix (used if update_type='ekf')
%       h_jacob     Function handle to generator for jacobian of 
%                   measurement matrix (used if msmt_type='ekf')
%
% Outputs:
%       x_pred      N x 1 predicted target state  (x_{k+1|k} if current 
%                   time index is k)
%       x_est       N x 1 estimated target state  (x_{k|k} if current time
%                   index is k)
%       P_pred      N x N predicted target state covariance matrix
%       P_est       N x N estimated target state covariance matrix
%
% Nicholas O'Donoughue
% 9 Dec 2021

%% Parse Inputs

if nargin < 8 || isempty(update_type)
    update_type = 'kf';
end
if nargin < 9 || isempty(msmt_type)
    msmt_type = 'kf';
end
valid_types = {'kf','ekf'};
assert(ismember(lower(update_type),valid_types),'Invalid KF update type.');
assert(ismember(lower(msmt_type),valid_types),'Invalid KF msmt type.');

do_coast = isempty(z);
N = size(x_pred,1); % Number of states

assert(size(P_pred,1)==N && size(P_pred,2)==N,'Dimension mismatch between target state and error covariance.');
assert(strcmpi(update_type,'ekf')==1 || (size(F,1)==N && size(F,2)==N),'Dimension mismatch for state transition matrix.');
assert(size(Q,1)==N && size(Q,2)==N,'Dimension mismsatch between target and process noise covariance.');

if ~do_coast
    M = size(z,1);      % Number of measurements
    assert(size(R,1)==M && size(R,2)==M,'Dimension mismatch between measurement and measurement noise covariance.');
    assert(strcmpi(msmt_type,'ekf')==1 || (size(H,1)==M && size(H,2)==N),'Dimension mismatch for measurement matrix.');
end




%% State Estimation Step
if do_coast
    % When coasting the measurement, there is no Kalman update.
    % x_est is just x_pred, and P_est is P_pred.
    x_est = x_pred;
    P_est = P_pred;
else
    % Compute the Kalman Gain and do update the current state estimate
    % based on the provided measurement
    
    % Parse KF or EKF model type
    switch lower(msmt_type)
        case 'kf'
            % H is a matrix
            z_pred = H*x_pred;
        case 'ekf'
            % H is a function handle, use it directly for z_pred, then
            % use the jacobian to compute H
            z_pred = H(x_pred);
            H = h_jacob(x_pred);
    end
    
    % Innovation Step, aka Measurement residual -- difference between 
    % predicted and actual msmt
    y = z - z_pred;

    % Compute error covariance for residual
    S = H*P_pred*H' + R;

    % Optimal Kalman Gain
    K = P_pred*H'/S;

    % Updated Target State Estimate
    x_est = x_pred + K*y;

    % Updated State Estimate Covariance
    P_est = (eye(N)-K*H)*P_pred;
end

%% State Prediction Step

% Parse KF or EKF model type
switch lower(update_type)
    case 'kf'
        % F is a matrix
        x_pred = F*x_est;
    case 'ekf'
        % F is a function handle, use it directly for x_pred, then
        % use the jacobian to compute H
        x_pred = F(x_est);
        F = f_jacob(x_est);
end

% Propagate forward the error covariance, and then add process noise
% covariance
P_pred = F*P_est*F' + Q;