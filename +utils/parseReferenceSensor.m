function [test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, num_sensors)
%[test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, num_sensors)
%
% Accepts a reference index setting (either None, a scalar integer, or a 
% 2 x N array of sensor pairs), and returns matching vectors for test 
% and reference indices.
%
% INPUTS:
%   ref_idx         Reference index setting
%   num_sensors     Number of available sensors
%
% OUTPUTS:
%   test_idx_vec    Vector of test sensor indices
%   ref_idx_vec     Vector of reference sensor indices
%
% Nicholas O'Donoughue
% 26 May 2021


%% Default Behavior
if isempty(ref_idx)
    % Default behavior is to generate all possible sensor pairs
    test_idx_vec = 1:num_sensors-1;
    ref_idx_vec = num_sensors * ones(size(test_idx_vec));
elseif strcmpi(ref_idx,'full')
    % Do the full set of N(N-1)/2 pairs
    full_set = nchoosek(1:num_sensors, 2);
    test_idx_vec = full_set(:,2)';
    ref_idx_vec = full_set(:,1)';
elseif isscalar(ref_idx)
    % Scalar reference provided, use all others as test sensors
    test_idx_vec = setdiff(1:num_sensors, ref_idx);
    ref_idx_vec = ref_idx * ones(size(test_idx_vec));
else
    % Explicit sensor pairs provided, parse them
    test_idx_vec = ref_idx(1, :);
    ref_idx_vec = ref_idx(2, :);
end
        
if size(test_idx_vec) ~= size(ref_idx_vec)
    warning('utils/parseReferenceSensor.m generated unequal test and reference vectors.  Check for bugs.');
end

return