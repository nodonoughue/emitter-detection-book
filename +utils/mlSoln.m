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
    else
        zz = 0;
    end
else
    yy = 0;
    zz = 0;
end

% Evaluate the likelihood function at each coordinate in the search space
A = zeros(numel(xx),numel(yy),numel(zz));
for idx_x = 1:numel(xx)
    for idx_y = 1:numel(yy)
        for idx_z = 1:numel(zz)
            switch nDim
                case 1
                    this_x = xx(idx_x);
                case 2
                    this_x = [xx(idx_x);yy(idx_y)];
                case 3
                    this_x = [xx(idx_x);yy(idx_y);zz(idx_z)];
            end
            A(idx_x,idx_y,idx_z) = ell(this_x);
        end
    end
end

% Find the peak
[~,idx_pk] = max(A(:));
[idx_pk_x,idx_pk_y,idx_pk_z] = ind2sub(size(A),idx_pk);
switch nDim
    case 1
        x_est = xx(idx_pk_x);
        x_grid = xx;
    case 2
        x_est = [xx(idx_pk_x);yy(idx_pk_y)];
        x_grid = {xx,yy};
    case 3
        x_est = [xx(idx_pk_x);yy(idx_pk_y);zz(idx_pk_z)];
        x_grid = {xx,yy,zz};
end
            