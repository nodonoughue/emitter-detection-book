function msmt_model = makeMeasurementModel(x_aoa, x_tdoa, x_fdoa, v_fdoa, ...
                                            tdoa_ref_idx, fdoa_ref_idx, state_space, ...
                                            R, least_square_fun, crlb_fun)
% makeMeasurementModel  Build a complete measurement model struct for the
%                       high-level tracker from sensor geometry.
%
% msmt_model = makeMeasurementModel(x_aoa, x_tdoa, x_fdoa, v_fdoa, ...
%                                   tdoa_ref_idx, fdoa_ref_idx, state_space)
% msmt_model = makeMeasurementModel(..., R)
% msmt_model = makeMeasurementModel(..., R, least_square_fun)
% msmt_model = makeMeasurementModel(..., R, least_square_fun, crlb_fun)
%
% INPUTS (geometry — same as the former makeMeasurementModel)
%   x_aoa        AOA sensor positions (2/3 x n_aoa), or []
%   x_tdoa       TDOA sensor positions (2/3 x n_tdoa), or []
%   x_fdoa       FDOA sensor positions (2/3 x n_fdoa), or []
%   v_fdoa       FDOA sensor velocities (2/3 x n_fdoa), or []
%   tdoa_ref_idx TDOA reference sensor index, or []
%   fdoa_ref_idx FDOA reference sensor index, or []
%   state_space  State space struct from makeMotionModel
%
% INPUTS (packaging — same as the former makeMsmtModel)
%   R                Measurement noise covariance (num_msmt x num_msmt), or []
%   least_square_fun Optional @(zeta, x_init) -> x_pos used by the initiator.
%                    Default: []
%   crlb_fun         Optional @(x_pos) -> (n x n) position CRLB.
%                    Default: []
%
% OUTPUTS
%   msmt_model  Struct with fields:
%                 z_fun            @(state_struct)  expected measurement
%                 h_fun            @(state_struct)  H matrix
%                 z_fun_raw        @(x)             measurement (raw vector)
%                 h_fun_raw        @(x)             H matrix (raw vector)
%                 R                measurement noise covariance (or [])
%                 state_space      state space struct
%                 least_square_fun function handle or []
%                 crlb_fun         function handle or []
%
% Nicholas O'Donoughue
% June 2025

if nargin < 8,  R                = []; end
if nargin < 9,  least_square_fun = []; end
if nargin < 10, crlb_fun         = []; end

% Sample the position/velocity components of the target state
x_pos = @(x) x(state_space.pos_idx, :);
if state_space.has_vel
    x_vel = @(x) x(state_space.vel_idx, :);
else
    x_vel = @(x) [];
end

if isempty(v_fdoa)
    v_fdoa = zeros(numel(state_space.vel_idx), 1);
end

%% Raw measurement function (accepts raw state vector x)
z_fun_raw = @(x) hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa - x_vel(x), ...
                                     x_pos(x), tdoa_ref_idx, fdoa_ref_idx);

%% Raw Jacobian function (accepts raw state vector x)
function h = buildH(x)
    [J, Jv] = hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_pos(x), ...
                               tdoa_ref_idx, fdoa_ref_idx, x_vel(x));
    [~, num_msmt, num_src] = size(J);
    h = zeros(num_msmt, state_space.num_states, num_src);
    h(:, state_space.pos_idx, :) = J.';
    if ~isempty(Jv)
        h(:, state_space.vel_idx, :) = Jv.';
    end
end

h_fun_raw = @buildH;

%% State-struct-aware wrappers (accept a State struct s)
z_fun = @(s) z_fun_raw(s.state);
h_fun = @(s) h_fun_raw(s.state);

msmt_model = struct('z_fun',            z_fun, ...
                    'h_fun',            h_fun, ...
                    'z_fun_raw',        z_fun_raw, ...
                    'h_fun_raw',        h_fun_raw, ...
                    'R',                R, ...
                    'state_space',      state_space, ...
                    'least_square_fun', least_square_fun, ...
                    'crlb_fun',         crlb_fun);
end
