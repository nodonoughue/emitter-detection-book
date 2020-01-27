function fd = dop(x1,v1,x2,v2,f)
% fd = dop(x1,v1,x2,v2,f)
%
% Given source and sensor at position x1 and x2 with velocity v1 and v2,
% compute the Doppler velocity shift
%
% INPUTS
%   x1      Position vector of N sources (nDim x N), in m
%   v1      Velocity vector of N sources (nDim x N), in m/s
%   x2      Position vector of M sensors (nDim x M), in m
%   v2      Velocity vector of M sensors (nDim x M), in m/s
%   f       Carrier frequency, in Hertz
%
% OUTPUTS
%   fd      Doppler shift for each source, sensor pair (N x M), in Hertz
%
% Nicholas O'Donoughue
% 1 July 2019

% Reshape inputs
[nDim,N] = size(x1);
[nDim2,M] = size(x2);

x1 = reshape(x1',N,1,nDim);
v1 = reshape(v1',[],1,nDim);
x2 = reshape(x2',1,M,nDim2);
v2 = reshape(v2',1,[],nDim2);

if nDim~=nDim2
    fprintf('Error: input dimensions do not match.');
    fd = [];
    return
end

% Unit vector from x1 to x2
u12 = (x2-x1)./sqrt(sum(abs(x2-x1).^2,3));
u21 = -u12;

% x1 velocity towards x2
vv1 = sum(v1.*u12,3);
vv2 = sum(v2.*u21,3);

% Sum of combined velocity
v = vv1 + vv2;

% Convert to Doppler
c = 3e8;
fd = f .* (1+ v./c);