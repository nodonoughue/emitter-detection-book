function r = rng(x1,x2)
% r = rng(x1,x2)
%
% Computes the range between two N-dimensional position vectors, using
% the Euclidean (L2) norm.
%
% Inputs:
%
%   x1          NxM1 matrix of N-dimensional positions
%   x2          NxM2 matrix of N-dimensional positions
%
% Outputs:
%
%   r           M1xM2 matrix of pair-wise ranges
%
% Nicholas O'Donoughue
% 1 July 2019

% Find the input dimensions and test for compatibility
[N1,M1] = size(x1);
[N2,M2] = size(x2);
if N1~=N2
    fprintf('Error; first dimension of x1 and x2 must match.\n');
end

% Reshape the inputs
x1 = reshape(x1,[N1,M1]);
x2 = reshape(x2,[N1,1,M2]);

% Compute the Euclidean Norm
r = reshape(sqrt(sum(abs(bsxfun(@minus,x1,x2)).^2,1)),[M1,M2]);