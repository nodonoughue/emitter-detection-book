function ddop = dopDiff(x0,v0,x1,v1,x2,v2,f)
% ddop = dopDiff(x0,v0,x1,v1,x2,v2,f)
%
% Computes the difference in Doppler shift between reference and test 
% sensor pairs and a source.  The two sets of sensor positions (x1 and x2)
% must have the same size.  Corresponding pairs will be compared.
%
% INPUTS:
%   x0      Position vector for N sources (nDim x N), in m
%   v0      Velocity vector for N sources (nDim x N), in m/s
%   x1      Position vector for M reference sensors (nDim x M), in m
%   v1      Velocity vector for M reference sensors (nDim x M), in m/s
%   x2      Position vector for M test sensors (nDim x M), in m
%   v2      Velocity vector for M test sensors (nDim x M), in m/s
%   f       Carrier frequency, in Hertz
%
% OUTPUTS:
%   ddop    Different Doppler shift (N x M), in Hertz
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute Doppler velocity from reference to each set of test positions
dop1 = utils.dop(x0,v0,x1,v1,f); % N x M
dop2 = utils.dop(x0,v0,x2,v2,f); % N x M

% Doppler difference
ddop = dop2 - dop1;
