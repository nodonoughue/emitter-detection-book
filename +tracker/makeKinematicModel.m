function [F, Q, state_space] = makeKinematicModel(model_type, num_dims, process_covar)
%




%% Parse Inputs

% num_dims, if not specified, defaults to 3
if nargin < 2 || isempty(num_dims)
    num_dims = 3;
end

if nargin < 3 || isempty(process_covar)
    process_covar = eye(num_dims);
elseif isscalar(process_covar)
    process_covar = process_covar*eye(num_dims);
end

%% Initialize Output
state_space = struct('num_dims',num_dims,...
                     'num_states',[],...
                     'has_pos',true,...
                     'has_vel',[],...
                     'pos_idx',[],...
                     'vel_idx',[]);

%% Define Kinematic Model and Process Noise
switch lower(model_type)
    case {'cv', 'constant velocity'}
        % Position and Velocity are tracked states
        % Acceleration is assumed zero-mean Gaussian
        state_space.num_states = 2*num_dims;
        state_space.pos_idx = 1:num_dims;
        state_space.vel_idx = num_dims + (1:num_dims);
        state_space.has_vel = true;

        F = @(t) [eye(num_dims), t*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims)];

        Q = @(t) [.25*t^4*process_covar, .5*t^3*process_covar;
                  .5*t^3*process_covar,     t^2*process_covar]; 
    case {'ca', 'constant acceleration'}
        % Position, Velocity, and Acceleration are tracked states
        % Acceleration is assumed to have non-zero-mean Gaussian
        % distribution
        %
        % State model is:
        % [px, py, pz, vx, vy, vz, ax, ay, az]'
        state_space.num_states = 3*num_dims;
        state_space.pos_idx = 1:num_dims;
        state_space.vel_idx = num_dims + (1:num_dims);
        state_space.has_vel = true;

        F = @(t) [eye(num_dims), t*eye(num_dims), .5*t^2*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims), t*eye(num_dims);
                  zeros(num_dims, 2*num_dims), eye(num_dims)];

        % Process noise covariance
        % This assumes that accel_var is the same in all dimensions
        Q = @(t) [.25*t^4*process_covar, .5*t^3*process_covar, .5*t^2*process_covar;
                   .5*t^3*process_covar,    t^2*process_covar,      t*process_covar;
                   .5*t^2*process_covar,      t*process_covar,        process_covar]; 

    case {'cj', 'constant jerk'}
        % Position, Velocity, and Acceleration are tracked states
        % Acceleration is assumed to have non-zero-mean Gaussian
        % distribution
        %
        % State model is:
        % [px, py, pz, vx, vy, vz, ax, ay, az, jx, jy, jz]'
        state_space.num_states = 4*num_dims;
        state_space.pos_idx = 1:num_dims;
        state_space.vel_idx = num_dims + (1:num_dims);
        state_space.has_vel = true;

        F = @(t) [eye(num_dims), t*eye(num_dims), .5*t^2*eye(num_dims), 1/6*t^3*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims), t*eye(num_dims), .5*t^2*eye(num_dims);
                  zeros(num_dims, 2*num_dims), eye(num_dims), t*eye(num_dims);
                  zeros(num_dims, 3*num_dims), eye(num_dims)];

        % Process noise covariance
        % This assumes that accel_var is the same in all dimensions
        Q = @(t) [t^7/252*process_covar, t^6/72*process_covar, t^5/30*process_covar, t^4/24*process_covar;
                  t^6/72*process_covar,  t^5/20*process_covar, t^4/8*process_covar,  t^3/6*process_covar;
                  t^5/30*process_covar,  t^4/8*process_covar,  t^3/3*process_covar,  t^2/2*process_covar; 
                  t^4/24*process_covar,  t^3/6*process_covar,  t^2/2*process_covar,      t*process_covar];

    % Implementation of the aerodynamic and ballistic models is left to
    % readers as an exercise
    case {'brv','ballistic reentry'}
        error('%s kinematic model not yet implemented.', model_type);
    case {'marv', 'maneuvering reentry'}
        error('%s kinematic model not yet implemented.', model_type);
    case {'aero'}
        error('%s kinematic model not yet implemented.', model_type);
    case {'ballistic'}
        error('%s kinematic model not yet implemented.', model_type);
    otherwise
        error('%s kinematic model option not recognized.', model_type);
end
