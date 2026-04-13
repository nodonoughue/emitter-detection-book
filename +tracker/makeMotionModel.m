function model = makeMotionModel(model_type, num_dims, process_covar, varargin)
% makeMotionModel  Create a unified motion model struct for use with the
%                  high-level tracker functions (predictState, runTracker, etc.).
%
% model = makeMotionModel(model_type, num_dims, process_covar)
% model = makeMotionModel('ct', num_dims, process_covar, process_covar_omega)
%
% INPUTS
%   model_type         String: 'cv', 'ca', 'cj', 'ballistic', 'ct'
%   num_dims           Number of spatial dimensions (default: 3)
%   process_covar      Scalar or (num_dims x num_dims) process noise covariance
%   process_covar_omega (CT only) Turn-rate spectral density [rad^2/s^3] (default: 1.0)
%
% OUTPUTS
%   model   Struct with fields common to all model types:
%             state_space   – state_space struct
%             is_linear     – true for CV/CA/CJ/Ballistic; false for CT
%             q_fun         – @(dt) process noise matrix
%           Fields for linear models (is_linear = true):
%             f_fun         – @(dt) transition matrix
%           Fields for ballistic model:
%             b_fun         – @(dt) control input vector u = [0;0;½g·dt²; 0;0;g·dt]
%           Fields for nonlinear models (is_linear = false):
%             f_fun_ekf     – @(x, dt) nonlinear transition function
%             jacobian_fun  – @(x, dt) EKF Jacobian
%
% Nicholas O'Donoughue
% June 2025

if nargin < 2 || isempty(num_dims)
    num_dims = 3;
end

if nargin < 3 || isempty(process_covar)
    process_covar = 1;
end

switch lower(model_type)
    case {'ct', 'constant_turn', 'constant turn'}
        % CT: nonlinear model
        if ~isempty(varargin)
            process_covar_omega = varargin{1};
        else
            process_covar_omega = 1.0;
        end

        [f_ct, g_ct, q_fun, ss] = make_constant_turn_model(num_dims, process_covar, process_covar_omega);

        model = struct('state_space',    ss, ...
                       'is_linear',      false, ...
                       'q_fun',          q_fun, ...
                       'f_fun',          [], ...
                       'f_fun_ekf',      f_ct, ...
                       'jacobian_fun',   g_ct);

    case {'ctra', 'constant_turn_rate_acceleration', 'constant turn rate acceleration'}
        % CTRA: nonlinear model with explicit acceleration states
        if ~isempty(varargin)
            process_covar_omega = varargin{1};
        else
            process_covar_omega = 1.0;
        end

        [f_ctra, g_ctra, q_fun, ss] = make_ctra_model(num_dims, process_covar, process_covar_omega);

        model = struct('state_space',    ss, ...
                       'is_linear',      false, ...
                       'q_fun',          q_fun, ...
                       'f_fun',          [], ...
                       'f_fun_ekf',      f_ctra, ...
                       'jacobian_fun',   g_ctra);

    otherwise
        % CV, CA, CJ, Ballistic: linear (or gravity-bias linear) model
        [F, Q, ss] = make_kinematic_model(model_type, num_dims, process_covar);

        model = struct('state_space',   ss, ...
                       'is_linear',     true, ...
                       'f_fun',         F, ...
                       'q_fun',         Q);

        % Ballistic: attach a control-input function for the gravity bias term.
        % u(dt) = [0;0;½g·dt²; 0;0;g·dt] — applied to mean only, not covariance.
        if strcmpi(model_type, 'ballistic')
            g   = -9.80665;               % m/s^2, downward
            p3  = ss.pos_idx(end);        % z-position index (3)
            v3  = ss.vel_idx(end);        % z-velocity index (6)
            n_s = ss.num_states;          % 6
            e_p = zeros(n_s, 1); e_p(p3) = 1;
            e_v = zeros(n_s, 1); e_v(v3) = 1;
            model.b_fun = @(dt) e_p * (0.5*g*dt^2) + e_v * (g*dt);
        end
end

end % makeMotionModel


%% -------------------------------------------------------------------------
function [F, Q, state_space] = make_kinematic_model(model_type, num_dims, process_covar)
% make_kinematic_model  Build transition and process-noise function handles
%                       for linear kinematic motion models (CV, CA, CJ, Ballistic).

if nargin < 2 || isempty(num_dims)
    num_dims = 3;
end

if nargin < 3 || isempty(process_covar)
    process_covar = eye(num_dims);
elseif isscalar(process_covar)
    process_covar = process_covar * eye(num_dims);
end

state_space = struct('num_dims',  num_dims, ...
                     'num_states', [], ...
                     'has_pos',   true, ...
                     'has_vel',   [], ...
                     'pos_idx',   [], ...
                     'vel_idx',   []);

switch lower(model_type)
    case {'cv', 'constant velocity'}
        state_space.num_states = 2*num_dims;
        state_space.pos_idx    = 1:num_dims;
        state_space.vel_idx    = num_dims + (1:num_dims);
        state_space.has_vel    = true;

        F = @(t) [eye(num_dims), t*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims)];

        Q = @(t) [.25*t^4*process_covar, .5*t^3*process_covar;
                  .5*t^3*process_covar,     t^2*process_covar];

    case {'ca', 'constant acceleration'}
        state_space.num_states = 3*num_dims;
        state_space.pos_idx    = 1:num_dims;
        state_space.vel_idx    = num_dims + (1:num_dims);
        state_space.accel_idx  = 2*num_dims + (1:num_dims);
        state_space.has_vel    = true;

        F = @(t) [eye(num_dims), t*eye(num_dims), .5*t^2*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims), t*eye(num_dims);
                  zeros(num_dims, 2*num_dims), eye(num_dims)];

        Q = @(t) [.25*t^4*process_covar, .5*t^3*process_covar, .5*t^2*process_covar;
                   .5*t^3*process_covar,    t^2*process_covar,      t*process_covar;
                   .5*t^2*process_covar,      t*process_covar,        process_covar];

    case {'cj', 'constant jerk'}
        state_space.num_states = 4*num_dims;
        state_space.pos_idx    = 1:num_dims;
        state_space.vel_idx    = num_dims + (1:num_dims);
        state_space.accel_idx  = 2*num_dims + (1:num_dims);
        state_space.jerk_idx   = 3*num_dims + (1:num_dims);
        state_space.has_vel    = true;

        F = @(t) [eye(num_dims),           t*eye(num_dims),    .5*t^2*eye(num_dims), 1/6*t^3*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims),      t*eye(num_dims),      .5*t^2*eye(num_dims);
                  zeros(num_dims,2*num_dims), eye(num_dims),    t*eye(num_dims);
                  zeros(num_dims,3*num_dims), eye(num_dims)];

        Q = @(t) [t^7/252*process_covar, t^6/72*process_covar, t^5/30*process_covar, t^4/24*process_covar;
                  t^6/72*process_covar,  t^5/20*process_covar, t^4/8*process_covar,  t^3/6*process_covar;
                  t^5/30*process_covar,  t^4/8*process_covar,  t^3/3*process_covar,  t^2/2*process_covar;
                  t^4/24*process_covar,  t^3/6*process_covar,  t^2/2*process_covar,      t*process_covar];

    case {'ballistic'}
        % Same linear CV dynamics as 'cv'.  Gravity enters as a deterministic
        % control input u(dt) attached as b_fun by the calling makeMotionModel.
        if num_dims ~= 3
            error('Ballistic model requires num_dims = 3 (z is the vertical axis).');
        end

        state_space.num_states = 2*num_dims;
        state_space.pos_idx    = 1:num_dims;
        state_space.vel_idx    = num_dims + (1:num_dims);
        state_space.has_vel    = true;

        F = @(t) [eye(num_dims), t*eye(num_dims);
                  zeros(num_dims,num_dims), eye(num_dims)];

        Q = @(t) [.25*t^4*process_covar, .5*t^3*process_covar;
                  .5*t^3*process_covar,     t^2*process_covar];

    case {'brv', 'ballistic reentry'}
        error('%s kinematic model not yet implemented.', model_type);
    case {'marv', 'maneuvering reentry'}
        error('%s kinematic model not yet implemented.', model_type);
    case {'aero'}
        error('%s kinematic model not yet implemented.', model_type);
    otherwise
        error('%s kinematic model option not recognized.', model_type);
end

end % make_kinematic_model


%% -------------------------------------------------------------------------
function [f_fun, g_fun, q_fun, state_space] = make_constant_turn_model(num_dims, process_covar, process_covar_omega)
% make_constant_turn_model  Build function handles for the Constant-Turn (CT) model.
%
% State layout:  [px, py, pz, vx, vy, vz, omega]   (7 states, num_dims must be 3)

if nargin < 1 || isempty(num_dims)
    num_dims = 3;
end
if num_dims ~= 3
    error('make_constant_turn_model: only num_dims = 3 (yaw-only CT) is supported.');
end

if nargin < 2 || isempty(process_covar)
    process_covar = eye(num_dims);
elseif isscalar(process_covar)
    process_covar = process_covar * eye(num_dims);
end

if nargin < 3 || isempty(process_covar_omega)
    process_covar_omega = 1.0;
end

num_states = 2*num_dims + 1;
state_space = struct('num_dims',   num_dims, ...
                     'num_states', num_states, ...
                     'has_pos',    true, ...
                     'has_vel',    true, ...
                     'is_linear',  false, ...
                     'pos_idx',    1:num_dims, ...
                     'vel_idx',    num_dims + (1:num_dims), ...
                     'omega_idx',  num_states);

f_fun = @(x, dt) ct_transition(x, dt);
g_fun = @(x, dt) ct_jacobian(x, dt);

q_cv  = @(t) [.25*t^4*process_covar, .5*t^3*process_covar;
               .5*t^3*process_covar,     t^2*process_covar];
q_fun = @(t) blkdiag(q_cv(t), process_covar_omega * t);

end % make_constant_turn_model


%% -------------------------------------------------------------------------
function x_next = ct_transition(x, dt)
% Exact discrete-time CT transition for state [px,py,pz,vx,vy,vz,omega].
% z-axis propagates as constant velocity (yaw-only rotation in x-y plane).

omega = x(7);
odt   = omega * dt;

if abs(odt) < 1e-6
    s_o  = dt;
    c1_o = 0.5 * odt * dt;
    c    = 1 - 0.5*odt^2;
    s    = odt;
else
    s_o  = sin(odt) / omega;
    c1_o = (1 - cos(odt)) / omega;
    c    = cos(odt);
    s    = sin(odt);
end

x_next = zeros(size(x));
x_next(1) = x(1) + x(4)*s_o  - x(5)*c1_o;
x_next(2) = x(2) + x(4)*c1_o + x(5)*s_o;
x_next(3) = x(3) + x(6)*dt;
x_next(4) = x(4)*c  - x(5)*s;
x_next(5) = x(4)*s  + x(5)*c;
x_next(6) = x(6);
x_next(7) = omega;

end % ct_transition


%% -------------------------------------------------------------------------
function G = ct_jacobian(x, dt)
% Analytical 7x7 Jacobian of ct_transition w.r.t. the state vector.

vx    = x(4);
vy    = x(5);
omega = x(7);
odt   = omega * dt;

if abs(odt) < 1e-6
    s_o      = dt;
    c1_o     = 0;
    c        = 1;
    s        = 0;
    ds_o_do  = 0;
    dc1_o_do = 0.5 * dt^2;
    dvx_do   = -vy * dt;
    dvy_do   =  vx * dt;
else
    s_o      = sin(odt) / omega;
    c1_o     = (1 - cos(odt)) / omega;
    c        = cos(odt);
    s        = sin(odt);
    ds_o_do  = dt*c/omega - s/omega^2;
    dc1_o_do = dt*s/omega - (1-c)/omega^2;
    dvx_do   = -vx*dt*s - vy*dt*c;
    dvy_do   =  vx*dt*c - vy*dt*s;
end

dpx_do = vx*ds_o_do  - vy*dc1_o_do;
dpy_do = vx*dc1_o_do + vy*ds_o_do;

G = zeros(7, 7);
G(1,1) = 1;  G(1,4) = s_o;   G(1,5) = -c1_o; G(1,7) = dpx_do;
G(2,2) = 1;  G(2,4) = c1_o;  G(2,5) = s_o;   G(2,7) = dpy_do;
G(3,3) = 1;  G(3,6) = dt;
G(4,4) = c;  G(4,5) = -s;    G(4,7) = dvx_do;
G(5,4) = s;  G(5,5) = c;     G(5,7) = dvy_do;
G(6,6) = 1;
G(7,7) = 1;

end % ct_jacobian


%% -------------------------------------------------------------------------
function [f_fun, g_fun, q_fun, state_space] = make_ctra_model(num_dims, process_covar, process_covar_omega)
% make_ctra_model  Build function handles for the Constant-Turn-Rate-and-
%                  Acceleration (CTRA) model.
%
% State layout:
%   2D: [px, py, vx, vy, ax, ay, omega]              — 7 states
%   3D: [px, py, pz, vx, vy, vz, ax, ay, az, omega]  — 10 states
%
% The turn rate omega (rad/s) and world-frame acceleration (ax, ay[, az])
% are tracked states.  Velocity rotates in the x-y plane at rate omega;
% the acceleration vector acts as an additive forcing term on position and
% velocity.  In 3D the z-axis propagates as constant-acceleration (CA).

if nargin < 1 || isempty(num_dims)
    num_dims = 3;
end
if num_dims ~= 2 && num_dims ~= 3
    error('make_ctra_model: num_dims must be 2 or 3.');
end
if nargin < 2 || isempty(process_covar)
    process_covar = eye(num_dims);
elseif isscalar(process_covar)
    process_covar = process_covar * eye(num_dims);
end
if nargin < 3 || isempty(process_covar_omega)
    process_covar_omega = 1.0;
end

num_states = 3*num_dims + 1;
state_space = struct('num_dims',   num_dims, ...
                     'num_states', num_states, ...
                     'has_pos',    true, ...
                     'has_vel',    true, ...
                     'is_linear',  false, ...
                     'pos_idx',    1:num_dims, ...
                     'vel_idx',    num_dims   + (1:num_dims), ...
                     'accel_idx',  2*num_dims + (1:num_dims), ...
                     'omega_idx',  num_states);

f_fun = @(x, dt) ctra_transition(x, dt, num_dims);
g_fun = @(x, dt) ctra_jacobian(x, dt, num_dims);

q_kin = @(t) [.25*t^4*process_covar, .5*t^3*process_covar, .5*t^2*process_covar;
               .5*t^3*process_covar,    t^2*process_covar,      t*process_covar;
               .5*t^2*process_covar,      t*process_covar,        process_covar];
q_fun = @(t) blkdiag(q_kin(t), process_covar_omega * t);

end % make_ctra_model


%% -------------------------------------------------------------------------
function x_next = ctra_transition(x, dt, n)
% Exact discrete-time CTRA transition.
% State: [pos(n) | vel(n) | accel(n) | omega]

omega = x(end);
odt   = omega * dt;

if abs(odt) < 1e-6
    sow = dt;
    com = 0.5 * odt * dt;   % (1-cos(odt))/omega ≈ omega*dt²/2
    c   = 1 - 0.5*odt^2;
    s   = odt;
else
    sow = sin(odt) / omega;
    com = (1 - cos(odt)) / omega;
    c   = cos(odt);
    s   = sin(odt);
end

vx = x(n+1);   vy = x(n+2);
ax = x(2*n+1); ay = x(2*n+2);

x_next = x;   % accel and omega are unchanged

x_next(1) = x(1) + sow*vx - com*vy + 0.5*dt^2*ax;   % px'
x_next(2) = x(2) + com*vx + sow*vy + 0.5*dt^2*ay;   % py'
x_next(n+1) =  c*vx - s*vy + dt*ax;                  % vx'
x_next(n+2) =  s*vx + c*vy + dt*ay;                  % vy'

if n == 3
    vz = x(n+3);
    az = x(2*n+3);
    x_next(3)   = x(3) + dt*vz + 0.5*dt^2*az;        % pz'
    x_next(n+3) = vz + dt*az;                         % vz'
end

end % ctra_transition


%% -------------------------------------------------------------------------
function G = ctra_jacobian(x, dt, n)
% Analytical Jacobian of ctra_transition w.r.t. the state vector.

ns    = 3*n + 1;
vx    = x(n+1);  vy = x(n+2);
omega = x(end);
odt   = omega * dt;

if abs(odt) < 1e-6
    sow      = dt;
    com      = 0.0;
    c        = 1;
    s        = 0;
    d_sow_do = 0;           % lim_{omega->0} d(sin(odt)/omega)/d(omega) = 0
    d_com_do = 0.5 * dt^2;  % lim_{omega->0} d((1-cos(odt))/omega)/d(omega) = dt²/2
else
    sow      = sin(odt) / omega;
    com      = (1 - cos(odt)) / omega;
    c        = cos(odt);
    s        = sin(odt);
    d_sow_do = dt*c/omega - s/omega^2;
    d_com_do = dt*s/omega - (1-c)/omega^2;
end

G = eye(ns);

% px row
G(1, n+1)   = sow;
G(1, n+2)   = -com;
G(1, 2*n+1) = 0.5*dt^2;
G(1, end)   = d_sow_do*vx - d_com_do*vy;

% py row
G(2, n+1)   = com;
G(2, n+2)   = sow;
G(2, 2*n+2) = 0.5*dt^2;
G(2, end)   = d_com_do*vx + d_sow_do*vy;

% vx row
G(n+1, n+1)   = c;
G(n+1, n+2)   = -s;
G(n+1, 2*n+1) = dt;
G(n+1, end)   = -dt*s*vx - dt*c*vy;

% vy row
G(n+2, n+1)   = s;
G(n+2, n+2)   = c;
G(n+2, 2*n+2) = dt;
G(n+2, end)   = dt*c*vx - dt*s*vy;

% 3D: z sub-system (pz', vz' rows)
if n == 3
    G(3,   n+3)   = dt;
    G(3,   2*n+3) = 0.5*dt^2;
    G(n+3, 2*n+3) = dt;
end

end % ctra_jacobian
