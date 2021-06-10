function cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
% cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
%
% Resample a covariance matrix based on a set of reference and test 
% indices.  This assumes a linear combination of the test and reference 
% vectors.  The output is an n_pair x n_pair covariance matrix for the 
% n_pair linear combinations.
%
% The measurements can be optionally weighted, to over-emphasize some
% measurements and de-emphasize others.
%
% In the resampled covariance matrix, the i,j-th entry is given
%    [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
%       where:  ai, aj are the i-th and j-th reference indices
%               bi, bj are the i-th and j-th test indices
%               C is the input covariance matrix
%
% Any indices that are NaN will be ignored, to represent a single-sensor
% measurement, such as AoA (which does not need a reference sensor
% measurement against which to compare), or a noise-free measurement with
% no error.
%
% If the third input, ref_idx, is missing or empty, then the second input,
% test_idx, will be passed to utils.parseReferenceSensor to generate
% matching test and reference vectors.
%
% INPUTS:
%   cov         NxN covariance matrix of individual sensor measurements
%   test_idx    n_pair x 1 vector of test sensor indices
%   ref_idx     n_pair x 1 vector of reference sensor indices [Optional]
%   test_wts    Optional n_pair x 1 vector of test measurement weights
%   ref_wts     Optional n_pair x 1 vector of reference measurement weights
%
% OUTPUTS:
%   cov_out     n_pair x n_pair output covariance matrix of sensor
%               measurement pairs
%
% Nicholas O'Donoughue
% 25 May 2021

%% Input handling
% Parse array sizes and indices
n_sensor = size(cov, 1);

% Handle test/reference inputs
if nargin < 3 || isempty(ref_idx)
    [test_idx, ref_idx] = utils.parseReferenceSensor(test_idx, n_sensor);
end

% Parse output size
n_test = numel(test_idx);
n_ref = numel(ref_idx);
n_out = max(n_test, n_ref);

if n_test > 1 && n_ref > 1 && n_test ~= n_ref
    error(strcat("Error calling covariance matrix resample. Reference and", ...
                 " test vectors must have the same shape."))
end

if any(test_idx > n_sensor) || any(ref_idx > n_sensor)
	error(strcat("Error calling covariance matrix resample. Indices exceed", ...
                 " the dimensions of the covariance matrix."))
end

% Parse sensor weights
do_test_wt = ~(nargin < 4 || isempty(test_wts));
do_ref_wt = ~(nargin < 5 || isempty(ref_wts));

if do_test_wt
    n_test_wt = numel(test_wts);
end

if do_ref_wt
    n_ref_wt = numel(ref_wts);
end

% Initialize output
cov_out = zeros(n_out, n_out);

a_i_wt = 1;
a_j_wt = 1;
b_i_wt = 1;
b_j_wt = 1;

% Step through reference sensors
for idx_row = 1:n_out
    % Parse sensor indices.  The mod commands seamlessly handle scalar
    % inputs
    a_i = test_idx(1+mod(idx_row-1, n_test));
    b_i = ref_idx(1+mod(idx_row-1, n_ref));

    % Parse sensor weights
    if do_test_wt
        a_i_wt = test_wts(1+mod(idx_row-1, n_test_wt));
    end

    if do_ref_wt
        b_i_wt = ref_wts(1+mod(idx_row-1, n_ref_wt));
    end

    for idx_col =1:n_out
        % Parse sensor indices.  The mod commands seamlessly handle scalar
        % inputs.
        a_j = test_idx(1+mod(idx_col-1, n_test));
        b_j = ref_idx(1+mod(idx_col-1, n_ref));

        if do_test_wt
            a_j_wt = test_wts(1+mod(idx_col-1, n_test_wt));
        end
    
        if do_ref_wt
            b_j_wt = ref_wts(1+mod(idx_col-1, n_ref_wt));
        end
        
        % Parse Input covariances
        if isnan(b_i) || isnan(b_j)
            cov_bibj = 0;
        else
            cov_bibj = cov(b_i, b_j);
        end
        if isnan(a_i) || isnan(a_j)
            cov_aiaj = 0;
        else
            cov_aiaj = cov(a_i, a_j);
        end
        if isnan(a_i) || isnan(b_j)
            cov_aibj = 0;
        else
            cov_aibj = cov(a_i, b_j);
        end
        if isnan(b_i) || isnan(a_j)
            cov_biaj = 0;
        else
            cov_biaj = cov(b_i, a_j);
        end
        
        %  [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
        cov_out(idx_row, idx_col) = b_i_wt * b_j_wt * cov_bibj + ...
                                    a_i_wt * a_j_wt * cov_aiaj - ...
                                    a_i_wt * b_j_wt * cov_aibj - ...
                                    b_i_wt * a_j_wt * cov_biaj;

    end
end