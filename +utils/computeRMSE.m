function rmse = computeRMSE(cov)
% rmse = computeRMSE(cov)
%
% Computes the root mean squared error (RMSE), defined for a 2D covariance
% matrix as the square root of the trace.
%
% This operates over the first two dimensions of C, and applies in parallel
% over any additional dimensions.
%
% Inputs:
%
% C         n_dim x n_dim covariance matrix (additional dimensions are
%           assumed to correspond to independent cases, and are computed
%           in turn)
%
% Outputs:
%
% rmse      root mean squared error
%
% Nicholas O'Donoughue
% 23 February 2022

% Check for multiple entries
c_dims = size(cov);
if numel(c_dims) > 2
    out_dims = c_dims(3:end);
    if numel(out_dims) == 1
        out_dims = [out_dims, 1];
    end
    n_cases = prod(out_dims);
else
    out_dims = [1, 1];
    n_cases = 1;
end

%% Compute RMSE
rmse = arrayfun(@(idx) sqrt(trace(cov(:,:,idx))), 1:n_cases);
rmse = reshape(rmse, out_dims);
