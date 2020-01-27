function af = compute_array_factor_ula(d_lam,N,psi,psi_0,el_pattern)
% af = compute_array_factor_ula(d_lam,N,psi,psi_0,el_pattern)
%
% Computes the array factor for a uniform linear array with specified 
% parameters.
%
% Inputs:
%
%   d_lam           Inter-element spacing (in wavelengths)
%   N               Number of array elements
%   psi             Incoming signal angle [radians]
%   psi_0           Steering angle [radians]
%   el_pattern      Optional element pattern (function handle that accepts
%                   psi and returns the individual element amplitude).
%
% Outputs:
%   af              Function handle that accepts a steering angle psi_0
%                   and test angle psi, both defined in radians, and
%                   returns 
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 4 || isempty(psi_0)
    psi_0 = pi/2; % Broadside is pi/2
end

if nargin < 5 || isempty(el_pattern)
    el_pattern = @(psi) 1;
end

% Build the Array Pattern
af = abs(sin(N*pi*d_lam*(sin(psi)-sin(psi_0)))...
    ./(N*(sin(pi*d_lam*(sin(psi)-sin(psi_0))))));

% Look for grating lobes
mask = abs(mod(d_lam*(sin(psi)-sin(psi_0))+.5,1)-.5)<1e-6;
af(mask) = 1;

% Apply the element pattern
el = el_pattern(psi);
af = af .* el;

