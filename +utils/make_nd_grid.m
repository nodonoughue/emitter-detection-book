function [x_set, dims] = make_nd_grid(ctr,max_offset,spacing)
% x_set = make_nd_grid(ctr,max_offset,spacing)
%
% Accepts a center value, maximum offset, and spacing vectors, of arbitrary
% length n_dim, and generates a search grid over n_dim dimensions.
%
% The returned matrix x_set has n_dim columns and N rows, where N is the
% product of the search vector in each dimension.
%    n_elements = 1 + 2*max_offset./spacing.
%    N = prod(n_elements)
%
% Inputs ctr, max_offset, and spacing must either be scalar, or have a
% common number of elements.  They are assumed to be vectors (shape is not
% preserved).
%
% Includes some error checking on the total array size allowable; which is
% currently set at no more than 10% of MATLAB's maximum array size, to be
% conservative.
%
% INPUTS:
%   ctr         Scalar or n_dim vector of center point for each dimension
%   max_offset  Scalar or n_dim vector of maximum search space offset from
%               ctr for each dimension
%   epsilon     Scalar or n_dim vector of search space resolution for each
%               dimension
%
% OUTPUTS:
%   x_set       n_dim x n_pt matrix of coordinates for n_pt test points
%   x_grid      n_dim x 1 cell array, each bearing the principal vector for
%               one dimension (to be used in plotting commands)
%
% Nicholas O'Donoughue
% 7 Nov 2021

%% Parse Inputs and Check for Max Array Size
n_dims = numel(ctr);

if numel(max_offset)==1
    max_offset = max_offset*ones(n_dims,1);
end

if numel(spacing)==1
    spacing = spacing * ones(n_dims,1);
end

assert(n_dims == numel(max_offset) && ...
       n_dims == numel(spacing),...
       'Search space dimensions do not match across specification of the center, search_size, and epsilon.');

n_elements = 1+ 2*max_offset(:)./spacing(:);

% Check Search Size
[user, ~] = memory;
maxElements = user.MaxPossibleArrayBytes/10; % Set a conservative limit,
                                             % since we'll need at least 2
assert(prod(n_elements) < maxElements, 'Search size is too large; MATLAB is likely to crash or become unresponsive.');

%% Initialize search space

% dims is a cell array of dimensions, each of which contains a vector of
% grid points along that dimension
dims = arrayfun(@(x,x_mx,n) x + x_mx * linspace(-1,1,n),  ctr(:), max_offset(:), n_elements(:), 'UniformOutput',false);

% Use meshgrid expansion; each element of dims_full is now a full n_dim 
% dimensioned grid for one of the elements of x
[x_grid{1:numel(dims)}] = ndgrid(dims{:}); %use comma-separated list expansion on both sides

% Rearrange to an n_dim x N matrix
x_grid_vec = cellfun(@(A) reshape(A,1,[]), x_grid, 'UniformOutput',false);
x_set = cell2mat(x_grid_vec(:)); % n_dim x N