function [x_est,A, x_grid] = mlSoln(ell,x_ctr,search_size,epsilon)
% function [x_est,A,x_grid] = mlSoln(ell,x_ctr,search_size,epsilon)
%
% Execute ML estimation through brute force computational methods.
%
% INPUTS:
%   ell          Function handle for the likelihood of a given position
%                must accept x_ctr (and similar sized vectors) as the sole
%                input.
%   x_ctr        Center position for search space (x, x/y, or z/y/z).
%   search_size  Search space size (same units as x_ctr)
%   epsilon      Search space resolution (same units as x_ctr)
%
% OUTPUTS:
%   x_est        Estimated minimum
%   A            Likelihood computed at each x position in the search space
%   x_grid       Set of x positions for the entire search space (M x N) for
%                N=1, 2, or 3.
%
% Nicholas O'Donoughue
% 1 July 2019

[x_set, x_grid] = utils.make_nd_grid(x_ctr, search_size, epsilon);

% rearrange to a matrix, where each column is
% Evaluate the likelihood function at each coordinate in the search space
A = ell(x_set);

% Find the peak
[~,idx_pk] = max(A(:));
x_est = x_set(:,idx_pk);
            