function z = measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, fdoa_ref_idx)
% z = measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, 
%                                                          fdoa_ref_idx)
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

% Construct component measurements
z_a =triang.measurement(x_aoa, x_source);
z_t =tdoa.measurement(x_tdoa, x_source,tdoa_ref_idx);
z_f =fdoa.measurement(x_fdoa, v_fdoa, x_source,fdoa_ref_idx);

% Combine into a single data vector
z = [z_a; z_t; z_f];