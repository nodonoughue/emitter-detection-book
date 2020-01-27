function dr = rngDiff(x0,x1,x2)
% dr = rngDiff(x0,x1,x2)
%
% Computes the difference in range from the reference input (x0) to each
% of the input vectors x1 and x2.  Difference is taken as the range from 
% x0 to each column in x2 minus the range from x0 to each column in x1, 
% pair-wise.
%
% The first dimension of all three inputs must match.  The dimensions of
% x1 and x2 determine the dimensions of the output dr.
%
% Inputs:
%
%   x0      Nx1 reference position
%   x1      NxM1 vector of test positions
%   x2      NxM2 vector of test positions
%
% Outputs:
%
%   dr      M1xM2 matrix of range differences
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute the range from the reference position to each of the set of
% test positions
r1 = utils.rng(x0,x1); % 1xM1
r2 = utils.rng(x0,x2); % 1xM2

% Take the difference, with appropriate dimension reshaping
dr = bsxfun(@minus,r2(:)',r1(:));