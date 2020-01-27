function af = compute_array_factor(v_fun,h,psi)
% af = compute_array_factor(v_fun,h,psi)
%
% Computes the array factor given a beamformer h, and array
% steering vector v (function handle) evaluated at a set of
% angles psi (in radians).
%
% The return is the response of the specified beamformer (h) to
% a plane wave (defined by v_fun) at various possible source angles.
% The outputs is in linear units of amplitude.
%
% Inputs:
%
%   v_fun       Function handle that returns N-element vector of complex
%               values
%   h           Beamformer fector (length N); will be normalized to peak
%               amplitude of 1
%   psi         Angles (in radians) over which to evaluate v_fun
%
% Outputs:
%
%   af          Linear array response of length N
%
% Nicholas O'Donoughue
% 1 July 2019

% Normalize the beamformer and make it a column vector
h = reshape(h,[],1)/max(abs(h(:)));

% Generate the steering vectors
vv = v_fun(psi); % should be numel(psi) x N

% Compute the inner product
af = vv'*h;


