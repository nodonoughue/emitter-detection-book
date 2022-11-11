function [zeta_mean, zeta_mean_full] = computeSampleMean(zeta)
% [zeta_mean, zeta_mean_full] = computeSampleMean(zeta)
%
% Computes the sample mean across the second dimension of zeta.
%
% Any additional dimensions (beyond the second) are assumed to be parallel
% trials, and are ignored.
%
% The first dimension is assumed to be the sample vector.
%
% INPUTS:
%   zeta        n_msmt x n_sample x [] matrix of measurements
%
% OUTPUTS:
%   zeta_mean   Sample mean across the second dimension
%   zeta_mean_full  Iterative sample mean (one for each sample)
%
% Nicholas O'Donoughue
% 28 February 2022

%% Parse Input
sum = cumsum(zeta,2);
K = 1:size(zeta,2);

%% Compute Sample Mean
zeta_mean_full = sum./(1:K);
zeta_mean = zeta_mean_full(:,end,:)

%% Reshape Output
dims = size(zeta);
if numel(dims)>2
    out_dims = [dims(1), dims(3:end)];
else
    out_dims = [dims(1), 1];
end

zeta_mean = reshape(zeta_mean,out_dims);