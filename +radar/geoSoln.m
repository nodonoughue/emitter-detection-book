function x = geoSoln(x_rdr,rho)
% x = geoSoln(x_rdr,rho)
%
% Uses geometry to estimate target position in cartesian coordinates. Does
% not consider range-rate measurements, and does not estimate target
% velocity.
%
% If there are four or more measurements, it assumes that they are:
%   range, range-rate, azimuth, and elevation, followed by extraneous
%   measurements (none after elevation are used)
%
% If there are three measurements, it assumes that they are:
%   range, azimuth, and elevation
%
% If there are fewer than three measurements, it throws an error.
%
% Inputs:
%   x_rdr               Radar positions [m]
%   rho                 Measurements [m]
% 
% Outputs:
%   x               Estimated source position
%
% Nicholas O'Donoughue
% 9 November 2022

% Parse inputs
[num_dim, ~] = size(x_rdr);
[num_msmt, ~] = size(rho);

assert(num_msmt > 2,'Not enough measurements to solve for target position.');

if num_msmt == 3
    range = rho(1,:);
    az = rho(2,:);
    el = rho(3,:);
else
    range = rho(1,:);
    az = rho(3,:);
    el = rho(4,:);
end

[e, n, u] = utils.aer2enu(az, el, range, 'rad', 'm');
x_tgt_off = [e(:), n(:), u(:)]'; % target position offset (rel. to radar)

if num_dim > 3
    % There are additional dimensions; make them NaN
    x_tgt_off(4:num_dim, :) = nan;
end

% Account for radar position
x = x_rdr + x_tgt_off;