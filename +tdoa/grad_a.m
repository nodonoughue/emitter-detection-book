function J = grad_a(x_tdoa, x_source, tdoa_ref_idx, ~)
% J = grad_a(x_tdoa, x_source,tdoa_ref_idx, alpha_tdoa)
%
% Returns the gradient of hybrid measurements, with sensor uncertainties,
% with respect to sensor measurement bias, alpha. Equation 6.28
%
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_tdoa          TDOA sensor positions
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%   do2DAOA
%   alpha_aoa
%   alpha_tdoa
%   alpha_fdoa
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%
% Nicholas O'Donoughue
% 31 October 2021

if nargin < 3 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx, size(x_tdoa,2));

%% Compute gradient
%  According to eq 6.29, the m-th row is c for every column for which the
%  m-th sensor is a test index, and -c for every column for which the m-th
%  sensor is a reference index.
num_tdoa = size(x_tdoa,2);
num_msmt = numel(test_idx_vec);
num_source = size(x_source,2);

% Convert idx vec from subscripts to indices
J = zeros(num_tdoa,num_msmt);
test_idx_set = sub2ind(size(J),test_idx_vec,1:num_msmt);
ref_idx_set = sub2ind(size(J),ref_idx_vec,1:num_msmt);

% Equation 6.29
J(test_idx_set) = 1;
J(ref_idx_set) = -1;

% Repeat for each source
J = repmat(J,1,1,num_source);