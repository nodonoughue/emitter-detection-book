function [zeta_update, K_update] = computeSampleMeanUpdate(zeta_old, zeta_new, K_original)
% [zeta_update, K_update] = computeSampleMeanUpdate(zeta_old, zeta_new, K_original)
%
% Computes the sample mean across the second dimension of zeta iteratively.
%
% INPUTS:
%   zeta_old        prior sample mean
%   zeta_new        new measurements to include in sample mean
%   K_original      number of samples in old sample mean
%
% OUTPUTS:
%   zeta_update     Updated sample mean
%   K_updated       Updated number of samples
%
% Nicholas O'Donoughue
% 28 February 2022

%% Parse Input
old_dims = size(zeta_old);
new_dims = size(zeta_new);
assert(old_dims(1)==new_dims(1),'Number of measurements must match between old and new sample mean.');

K_new = new_dims(2);
if numel(new_dims)>2
    new_dims_out = [new_dims(1),new_dims(3:end)];
else
    new_dims_out = [new_dims(1),1];
end

K_update = K_original + K_new;

%% Compute Sample Mean Update
zeta_innovation = reshape(sum(zeta_new,2)/K_update,new_dims_out);

zeta_update = zeta_original * K_original/K_update + zeta_innovation;