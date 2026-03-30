function msmt_model = makeMsmtModel(z_fun_raw, h_fun_raw, R, state_space, least_square_fun, crlb_fun)
% makeMsmtModel  Create a measurement model struct for the high-level tracker.
%
% msmt_model = makeMsmtModel(z_fun_raw, h_fun_raw, R, state_space)
% msmt_model = makeMsmtModel(z_fun_raw, h_fun_raw, R, state_space, least_square_fun)
% msmt_model = makeMsmtModel(z_fun_raw, h_fun_raw, R, state_space, least_square_fun, crlb_fun)
%
% The raw function handles z_fun_raw / h_fun_raw are those produced by
% tracker.makeMeasurementModel — they accept a raw state vector.
% This wrapper builds State-struct-aware versions for use with the tracker
% infrastructure (associateTracks, runTrackerStep, etc.).
%
% INPUTS
%   z_fun_raw        @(x) → measurement vector (x is num_states x 1)
%   h_fun_raw        @(x) → H matrix (num_msmt x num_states)
%   R                Measurement noise covariance matrix (num_msmt x num_msmt)
%   state_space      State space struct from makeMotionModel
%   least_square_fun Optional @(zeta, x_init) → [x_pos, ~] function used by the
%                    TwoPointInitiator to estimate a position from a measurement.
%                    Default: [] (initiator will not be able to convert measurements)
%   crlb_fun         Optional @(x_pos) → (num_dims x num_dims) position CRLB matrix.
%                    Used by initiateTracks to compute buffer/track covariances via
%                    a dedicated CRLB function (e.g. tdoa.computeCRLB) rather than
%                    the manual H'*R^{-1}*H FIM inversion, which can be numerically
%                    unstable for false-alarm positions near sensor singularities.
%                    Default: [] (fall back to manual FIM from h_fun)
%
% OUTPUTS
%   msmt_model  Struct with fields:
%                 z_fun            @(state_struct)  expected measurement vector
%                 h_fun            @(state_struct)  H matrix (num_msmt x num_states)
%                 R                measurement noise covariance
%                 state_space      reference to the state space
%                 least_square_fun function handle, or []
%                 crlb_fun         function handle, or []
%
% Nicholas O'Donoughue
% June 2025

if nargin < 5
    least_square_fun = [];
end

if nargin < 6
    crlb_fun = [];
end

% Wrap the raw function handles to accept State structs
z_fun = @(s) z_fun_raw(s.state);
h_fun = @(s) h_fun_raw(s.state);

msmt_model = struct('z_fun',            z_fun, ...
                    'h_fun',            h_fun, ...
                    'R',                R, ...
                    'state_space',      state_space, ...
                    'least_square_fun', least_square_fun, ...
                    'crlb_fun',         crlb_fun);
