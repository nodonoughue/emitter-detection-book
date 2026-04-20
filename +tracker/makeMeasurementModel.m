function msmt_model = makeMeasurementModel(x_aoa, x_tdoa, x_fdoa, v_fdoa, ...
                                            tdoa_ref_idx, fdoa_ref_idx, ...
                                            R, least_square_fun, crlb_fun)
% makeMeasurementModel  Build a complete measurement model struct for the
%                       high-level tracker from sensor geometry.
%
% msmt_model = makeMeasurementModel(x_aoa, x_tdoa, x_fdoa, v_fdoa, ...
%                                   tdoa_ref_idx, fdoa_ref_idx)
% msmt_model = makeMeasurementModel(..., R)
% msmt_model = makeMeasurementModel(..., R, least_square_fun)
% msmt_model = makeMeasurementModel(..., R, least_square_fun, crlb_fun)
%
% INPUTS (geometry)
%   x_aoa        AOA sensor positions (2/3 x n_aoa), or []
%   x_tdoa       TDOA sensor positions (2/3 x n_tdoa), or []
%   x_fdoa       FDOA sensor positions (2/3 x n_fdoa), or []
%   v_fdoa       FDOA sensor velocities (2/3 x n_fdoa), or []
%   tdoa_ref_idx TDOA reference sensor index, or []
%   fdoa_ref_idx FDOA reference sensor index, or []
%
% INPUTS (packaging)
%   R                Measurement noise covariance (num_msmt x num_msmt), or []
%   least_square_fun Optional @(zeta, x_init) -> x_pos used by the initiator.
%                    Default: []
%   crlb_fun         Optional @(x_pos) -> (n x n) position CRLB.
%                    Default: []
%
% OUTPUTS
%   msmt_model  Struct with fields:
%                 z_fun            @(state_struct)  expected measurement
%                 h_fun            @(state_struct)  H matrix (num_msmt x num_states)
%                 R                measurement noise covariance (or [])
%                 least_square_fun function handle or []
%                 crlb_fun         function handle or []
%
% The state_space is no longer stored in the struct; z_fun and h_fun extract
% pos_idx / vel_idx / num_states from the State struct passed at call time.
%
% Nicholas O'Donoughue
% June 2025

if nargin < 7,  R                = []; end
if nargin < 8,  least_square_fun = []; end
if nargin < 9,  crlb_fun         = []; end

msmt_model = struct('z_fun',            @compute_z, ...
                    'h_fun',            @compute_h, ...
                    'R',                R, ...
                    'least_square_fun', least_square_fun, ...
                    'crlb_fun',         crlb_fun);

    %% State-struct-aware measurement function
    function z = compute_z(s)
        ss        = s.state_space;
        x_pos_val = s.state(ss.pos_idx);
        if ss.has_vel
            x_vel_val = s.state(ss.vel_idx);
        else
            x_vel_val = [];
        end
        v_rel = v_fdoa;
        if isempty(v_rel) && ~isempty(x_fdoa)
            v_rel = zeros(size(x_vel_val));
        elseif ~isempty(v_rel)
            v_rel = v_rel - x_vel_val;
        end
        z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_rel, x_pos_val, ...
                                tdoa_ref_idx, fdoa_ref_idx);
    end

    %% State-struct-aware Jacobian function
    function h = compute_h(s)
        ss        = s.state_space;
        x_pos_val = s.state(ss.pos_idx);
        if ss.has_vel
            x_vel_val = s.state(ss.vel_idx);
        else
            x_vel_val = [];
        end
        [J, Jv] = hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_pos_val, ...
                                   tdoa_ref_idx, fdoa_ref_idx, x_vel_val);
        [~, num_msmt, num_src] = size(J);
        h = zeros(num_msmt, ss.num_states, num_src);
        h(:, ss.pos_idx, :) = J.';
        if ~isempty(Jv)
            h(:, ss.vel_idx, :) = Jv.';
        end
    end

end
