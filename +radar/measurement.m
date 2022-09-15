function rho = measurement(x_rdr, x_tgt, v_rdr, v_tgt, angle_dims)
% 
%
% Computed monostatic range and range rate (Doppler) measurements for a 
% set of bistatic transmitter/receiver pairs.
%
% If the velocity terms are all blank, then range only measurements are
% generated.  Otherwise range-rate measurements are generated, and returned
% stacked vertically after the range measurements.
%
% INPUTS:
%   x_rdr            nDim x nRdr array of radar coordinates
%   x_tgt            nDim x nTgt array of target coordinates
%   v_rdr            nDim x nRdr array of radar velocities
%   v_tgt            nDim x nTgt array of target velocities
%   angle_dims       Flag determining how many angle of arrival dimensions 
%                    to report (0 =  no angle data, 1 = azimuth, 2 =
%                    azimuth & elevation)
% OUTPUTS:
%   rho              nMsmt x nTgt matrix of range (and potentially
%                    range-rate and angle of arrival) measurements.
%                    Measurements are ordered range, range-rate, azimuth
%                    and elevation.
%
% Nicholas O'Donoughue
% 15 September 2022

do_doppler = nargin >= 2 && (~isempty(v_rdr) || ~isempty(v_tgt));

if nargin < 5 || ~exist('angle_dims','var') || isempty(angle_dims)
    angle_dims = 0;
end
assert(angle_dims >=0 && angle_dims <= 2,'Error parsing angle dimensions command, must be between 0 and 2.');

if do_doppler && isempty(v_rdr)
    v_rdr = zeros(size(x_rdr));
end

if do_doppler && isempty(v_tgt)
    v_tgt = zeros(size(x_tgt));
end

% Parse inputs
[nDim,nRdr] = size(x_rdr);
[~,nTgt] = size(x_tgt);

%% Bistatic Range
% Compute distance from each source position to each sensor
dx = reshape(x_tgt,nDim,1,nTgt) - reshape(x_rdr,nDim,nRdr);
rng = sqrt(sum(abs(dx).^2,1)); % 1 x nTx x nTgt

%% Bistatic Range Rate
if do_doppler
    % Compute range rate from range and velocity
    dv = reshape(v_tgt,nDim,1,nTgt) - reshape(v_rdr,nDim,nRdr);

    rng_rate = reshape(-sum(dv.*dx./rng,1),nRdr,nTgt);
else
    rng_rate = [];
end

% Append to Measurement Matrix
rho = cat(1,rng, rng_rate);

%% Receiver angle measurements
if angle_dims > 0
    % Check dimension output
    do2DAOA = (angle_dims == 2);
        
	% Generate measurements
    psi = triang.measurement(x_rx,x_tgt,do2DAOA);
    rho = cat(1,rho,psi);
end
