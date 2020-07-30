function [x_est,A,x_grid] = mlSoln(ell,x_ctr,search_size,epsilon)
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

nDim = numel(x_ctr);
if nDim <1 || nDim > 3
    error('Number of spatial dimensions must be between 1 and 3');
end

if numel(search_size)==1
    search_size = search_size*ones(nDim,1);
end

% Initialize search space
xx = x_ctr(1) + (-search_size(1):epsilon:search_size(1));
if nDim > 1
    yy = x_ctr(2) + (-search_size(2):epsilon:search_size(2));
    if nDim > 2
        zz = x_ctr(3) + (-search_size(3):epsilon:search_size(3));
        [XX,YY,ZZ] = ndgrid(xx,yy,zz);
        x_set = [XX(:),YY(:),ZZ(:)]';
        x_grid = {xx,yy,zz};
    else
        [XX,YY] = ndgrid(xx,yy);
        x_set = [XX(:),YY(:)]';
        x_grid = {xx,yy};
    end
else
    x_set = xx(:)';
    x_grid = xx;
end

% Evaluate the likelihood function at each coordinate in the search space
A = ell(x_set);

% Find the peak
[~,idx_pk] = max(A(:));
x_est = x_set(:,idx_pk);
            