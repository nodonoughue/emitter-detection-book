function m = makeMeasurement(time, zeta, sensor_info)
% makeMeasurement  Create a tracker measurement struct.
%
% m = makeMeasurement(time, zeta)
% m = makeMeasurement(time, zeta, sensor_info)
%
% INPUTS
%   time         Scalar timestamp [s]
%   zeta         Column vector of measured values
%   sensor_info  Optional struct containing sensor parameters (e.g. the
%                sensor position array x_sensor and covariance C used when
%                building a measurement model).  Default: []
%
% OUTPUTS
%   m   Struct with fields:
%         time         – scalar timestamp [s]
%         zeta         – (num_measurements x 1) measurement vector
%         sensor_info  – sensor parameter struct, or []
%
% Nicholas O'Donoughue
% June 2025

if nargin < 3
    sensor_info = [];
end

m = struct('time',        time, ...
           'zeta',        zeta(:), ...
           'sensor_info', sensor_info);
