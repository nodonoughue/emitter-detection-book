function z = measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, fdoa_ref_idxm, do2Daoa)
% z = measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, 
%                                                fdoa_ref_idx, do2Daoa)
%
% Computes hybrid measurements, for AOA, TDOA, and FDOA sensors.
%
% INPUTS:
%   x_aoa       nDim x nAOA array of sensor positions
%   x_tdoa      nDim x nTDOA array of TDOA sensor positions
%   x_fdoa      nDim x nFDOA array of FDOA sensor positions
%   v_fdoa      nDim x nFDOA array of FDOA sensor velocities
%   x_source    nDim x nSource array of source positions
%   tdoa_ref_idx    Index for reference TDOA sensor or 2 x nPair set of
%                   TDOA sensor pairing indices [optional]
%   fdoa_ref_idx    Index for reference FDOA sensor or 2 x nPair set of
%                   FDOA sensor pairing indices [optional]
%   do2Daoa     Boolean flag, if true then 2D AOA measurements will be
%               generated (azimuth and elevation)
%
% OUTPUTS:
%   z           nAoa + nTDOA + nFDOA - 2 x nSource array of measurements
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end
if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 8 || ~exist('do2Daoa','var')
    do2Daoa = true;
end

% Construct component measurements
if ~isempty(x_aoa)
    z_a =triang.measurement(x_aoa, x_source, do2Daoa);
else
    z_a = [];
end
if ~isempty(x_tdoa)
    z_t =tdoa.measurement(x_tdoa, x_source,tdoa_ref_idx);
else
z_t = [];
end
if ~isempty(x_fdoa)
    z_f =fdoa.measurement(x_fdoa, v_fdoa, x_source,fdoa_ref_idx);
else
    z_f = [];
end

% Combine into a single data vector
z = [z_a; z_t; z_f];