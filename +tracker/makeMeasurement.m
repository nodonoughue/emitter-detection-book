function m = makeMeasurement(msmt_model_or_time, state_or_zeta, time_or_msmt_model)
% makeMeasurement  Create a tracker measurement struct.
%
% Primary form — generate zeta from a state via the measurement model:
%   m = makeMeasurement(msmt_model, state, time)
%
% Explicit-zeta form — wrap a pre-computed measurement vector:
%   m = makeMeasurement(time, zeta)
%   m = makeMeasurement(time, zeta, msmt_model)
%
% INPUTS (primary form)
%   msmt_model   Measurement model struct from makeMeasurementModel
%   state        State struct from makeState (used as input to z_fun)
%   time         Scalar timestamp [s]
%
% INPUTS (explicit-zeta form)
%   time         Scalar timestamp [s]
%   zeta         Column vector of measured values
%   msmt_model   Optional measurement model struct (default: [])
%
% OUTPUTS
%   m   Struct with fields:
%         time       – scalar timestamp [s]
%         zeta       – (num_measurements x 1) measurement vector
%         msmt_model – measurement model struct, or []
%
% Nicholas O'Donoughue
% June 2025

if isstruct(msmt_model_or_time)
    % Primary form: (msmt_model, state, time)
    msmt_model = msmt_model_or_time;
    state      = state_or_zeta;
    time       = time_or_msmt_model;
    zeta       = msmt_model.z_fun(state);
else
    % Explicit-zeta form: (time, zeta) or (time, zeta, msmt_model)
    time  = msmt_model_or_time;
    zeta  = state_or_zeta;
    if nargin >= 3
        msmt_model = time_or_msmt_model;
    else
        msmt_model = [];
    end
end

m = struct('time',       time, ...
           'zeta',       zeta(:), ...
           'msmt_model', msmt_model);
