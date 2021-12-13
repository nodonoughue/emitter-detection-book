function [x_est,A,x_grid] = bestfix(pdfs,x_ctr,search_size,epsilon)
%
% INPUTS:
%   pdfs        Lx1 cell array of PDF functions, each of which represents
%               the probability that an input position (Ndim x 3 array)
%               is the true source position for one of the measurements.
%   x_ctr        Center position for search space (x, x/y, or z/y/z).
%   search_size  Search space size (same units as x_ctr)
%   epsilon      Search space resolution (same units as x_ctr)
%
% OUTPUTS:
%   x_est        Estimated position
%   A            Likelihood computed at each x position in the search space
%   x_grid       Set of x positions for the entire search space (M x N) for
%                N=1, 2, or 3.
%
%
% Based on the BESTFIX algorithm, invented in 1990 by Eric Hodson (R&D 
% Associates, now with Naval Postgraduate School).  Patent is believed to 
% be in the public domain.
%
% Ref:
%  Eric Hodson, "Method and arrangement for probabilistic determination of 
%  a target location," U.S. Patent US5045860A, 1990, 
%  https://patents.google.com/patent/US5045860A
%
%
% Nicholas A. O'Donoughue
% 29 July 2020

[x_set, x_grid] = utils.make_nd_grid(x_ctr, search_size, epsilon);

A = ones(1,size(x_set,2));

% Loop across the sensors
for idx = 1:numel(pdfs)
    A = A .* pdfs{idx}(x_set);
end

% Find the highest scoring position
[~,idx_pk] = max(A(:));
x_est = x_set(:,idx_pk);