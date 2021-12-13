function [x_aoa, x_tdoa, x_fdoa, v_fdoa, C_beta] = ...
               parsePosVelCovar(x_sensor, v_sensor, aoa_idx, tdoa_idx,...
                                fdoa_idx, C_pos, C_vel)
% [x_aoa, x_tdoa, x_fdoa, v_fdoa, C_beta] = ...
%               parsePosVelCovar(x_sensor, v_sensor, aoa_idx, tdoa_idx,...
%                                fdoa_idx, C_pos, C_vel)
%
% Parses a set of positions and velocities according to the provided
% indices, to generate separate sets of position for AOA sensors, position
% for TDOA sensors, and position and velocity for FDOA sensors.
%
% Optionally, provide a set of position and velocity covariance matrices
% for each sensor.  If provided, then the sensor indices are used to create
% a 2D covariance matrix for the aggregate parameter beta that collects
% sensor position and velocity across all sensor types:
%
%       beta0 = [xaoa; xtdoa; xfdoa; vfdoa]
%       Cbeta = covariance matrix of sensor pos/vel
%       beta ~ N(beta0,Cbeta)       
%
% To skip a sensor type (e.g. skip AOA for a TDOA/FDOA scenario), specify
% the indices for that sensor as an empty array.
%
% INPUTS:
%   x_sensor    n_dim x n_sensor array of sensor positions (or reported
%               positions)
%   v_sensor    n_dim x n_sensor array of sensor velocities (if no FDOA
%               sensors are needed, this input can be an empty array).
%   aoa_idx     n_aoa x 1 vector of indices for sensors that provide AOA
%               measurements (use an empty array to skip)
%   tdoa_idx    n_tdoa x 1 vector of indices for sensors that provide TDOA
%               measurements (use an empty array to skip)
%   fdoa_idx    n_fdoa x 1 vector of indices for sensors that provide FDOA
%               measurements (use an empty array to skip).
%   C_pos       n_dim x n_dim x n_sensor collection of position covariance
%               matrices for each sensor. If not provided, then no output
%               covariance matrix will be generated (an empty array will be
%               returned).
%   C_vel       n_dim x n_dim x n_sensor collection of position covariance
%               matrices for each sensor. Required if n_fdoa > 0 and C_pos
%               is provided.
%
% OUTPUTS:
%   x_aoa       n_dim x n_aoa array of AOA sensor positions
%   x_tdoa      n_dim x n_tdoa array of TDOA sensor positions
%   x_fdoa      n_dim x n_fdoa array of FDOA sensor positions
%   v_fdoa      n_dim x n_fdoa array of FDOA sensor velocities
%   C_beta      Square covariance matrix for sensor positions and
%               velocities, with n_dim * (n_aoa + n_tdoa + 2*n_fdoa) 
%               rows and columns.
%
% Nicholas O'Donoughue
% 5 Nov 2021

%% Parse Inputs

if nargin < 3 || ~exist('aoa_idx','var')
    aoa_idx = [];
end

if nargin < 4 || ~exist('tdoa_idx','var')
    tdoa_idx = [];
end

if nargin < 5 || ~exist('fdoa_idx','var')
    fdoa_idx = [];
end

assert(~isempty(aoa_idx) || ~isempty(tdoa_idx) || ~isempty(fdoa_idx),...
       'No sensors specified, please provide indices for at least one of AOA, TDOA, and FDOA.');

n_dim = size(x_sensor,1);
n_sensor = size(x_sensor,2);
n_vel = size(v_sensor,2);
n_aoa = numel(aoa_idx);
n_tdoa = numel(tdoa_idx);
n_fdoa = numel(fdoa_idx);

assert(~any(aoa_idx > n_sensor) && ~any(tdoa_idx > n_sensor) ...
       && ~any(fdoa_idx > n_sensor) && ~any(fdoa_idx > n_vel),...
       'One or more of the sensor indices exceeds the number of sensors provided.');

do_covar = nargin >= 6 && exist('C_pos','var') && ~isempty(C_pos);
do_fdoa_covar = nargin >= 7 && exist('C_vel','var') && ~isempty(C_vel);
if do_covar && n_fdoa > 0
    assert(do_fdoa_covar, 'Missing sensor velocity covariance.');
end

%% Parse Sensor Positions

if isempty(aoa_idx)
    x_aoa = [];
else
    x_aoa = x_sensor(:,aoa_idx);
end

if isempty(tdoa_idx)
    x_tdoa = [];
else
    x_tdoa = x_sensor(:,tdoa_idx);
end

if isempty(fdoa_idx)
    x_fdoa = [];
    v_fdoa = [];
else
    x_fdoa = x_sensor(:,fdoa_idx);
    v_fdoa = v_sensor(:,fdoa_idx);
end

%% Parse Sensor Covariances
if ~do_covar
    C_beta = [];
else
    C_beta = zeros(n_dim * (n_aoa + n_tdoa + 2*n_fdoa));

    % Fill AOA covariances
    for idx = 1:n_aoa
        n = aoa_idx(idx);               % sensor number
        this_C = squeeze(C_pos(:,:,n)); % sensor covariance

        % Assign to the idx-th block diagonal
        cov_idx = (idx-1)*n_dim + (1:n_dim);
        C_beta(cov_idx,cov_idx) = this_C;
    end

    % Fill TDOA covariances
    for idx = 1:n_tdoa
        n = tdoa_idx(idx);              % sensor number
        this_C = squeeze(C_pos(:,:,n)); % sensor covariance

        % Assign to the idx-th block diagonal
        cov_idx = (n_aoa+idx-1)*n_dim + (1:n_dim);
        C_beta(cov_idx,cov_idx) = this_C;

        % Check if this sensor is also an AOA sensor
        idx_aoa = find(aoa_idx == n);   % position of this sensor in the 
                                        % AOA set of sensors
        cov_idx_aoa = (idx_aoa-1)*n_dim + (1:n_dim);

        % Assign this_C to the cross terms for this sensor as reflected in
        % the AOA and TDOA sets
        C_beta(cov_idx,cov_idx_aoa) = this_C;
        C_beta(cov_idx_aoa,cov_idx) = this_C;
    end

    % Fill FDOA covariances
    for idx = 1:n_fdoa
        n = fdoa_idx(idx);                  % sensor number
        this_C = squeeze(C_pos(:,:,n));     % sensor covariance
        this_C_vel = squeeze(C_vel(:,:,n)); % sensor velocity covariance

        % Assign to the idx-th block diagonal
        cov_idx = (n_aoa+n_tdoa+idx-1)*n_dim + (1:n_dim);
        cov_idx_vel = cov_idx + n_fdoa;
        C_beta(cov_idx, cov_idx) = this_C;
        C_beta(cov_idx_vel, cov_idx_vel) = this_C_vel;

        % Check if this sensor is also an AOA sensor
        idx_aoa = find(aoa_idx == n);   % position of this sensor in the 
                                        % AOA set of sensors
        cov_idx_aoa = (idx_aoa-1)*n_dim + (1:n_dim);

        % Assign this_C to the cross terms for this sensor as reflected in
        % the AOA and TDOA sets
        C_beta(cov_idx,cov_idx_aoa) = this_C;
        C_beta(cov_idx_aoa,cov_idx) = this_C;

        % Check if this sensor is also a TDOA sensor
        idx_tdoa = find(tdoa_idx == n);   % position of this sensor in the 
                                          % TDOA set of sensors
        cov_idx_tdoa = (n_aoa+idx_tdoa-1)*n_dim + (1:n_dim);

        % Assign this_C to the cross terms for this sensor as reflected in
        % the AOA and TDOA sets
        C_beta(cov_idx,cov_idx_tdoa) = this_C;
        C_beta(cov_idx_tdoa,cov_idx) = this_C;
    end
end