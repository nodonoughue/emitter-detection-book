function cov_out = ensureInvertible(cov, epsilon)
% function cov_out = ensureInvertible(cov, epsilon)
%
% Check the input matrix for invertibility by finding the eigenvalues and
% checking that they are all >= a small value (epsilon).
%
% If any of the eigenvalues are too small, then a diagonal loading term
% is applied to ensure that the matrix is positive definite (all
% eigenvalues are >= epsilon).
%
% INPUTS:
%   cov     Input square covariance matrix.  If the input has >2
%           dimensions, then the process is repeated across the extra
%           dimensions.
%
%   epsilon (Optional) Specifies the minimum eigenvalue to use when
%           checking for invertibility. Default = 1e-10
%
% OUTPUTS:
%   cov_out Output covariance matrix, guaranteed to be invertible.
%
% Nicholas O'Donoughue
% 1 June 2021

% Check for epsilon input
if nargin < 2 || isempty(epsilon)
    epsilon = 1e-20;
end

% Check input dimensions
sz = size(cov);
assert(numel(sz) > 1, 'Input must have at least two dimensions.');
assert(sz(1) == sz(2), 'First two dimensions of input matrix must be equal.');
dim = sz(1);
if numel(sz) > 2
    n_matrices = prod(sz(3:end));
else
    n_matrices = 1;
end

cov_out = zeros(size(cov));
for idx_matrix = 1:n_matrices
   % Check min eigenvalue
   if min(eig(cov(:,:,idx_matrix))) < epsilon
      % We need to add a diagonal loading term, determine appropriate size
      d = epsilon;
      while min(eig(cov(:,:,idx_matrix) + d * eye(dim))) < epsilon
          d = d * 10;
      end
      
      % Apply diagonal loading
      cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix) + d * eye(dim);
   else
       % No diagonal loading
       cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix);
   end
end
