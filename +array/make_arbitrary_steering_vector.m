function [v,v_dot] = make_arbitrary_steering_vector(d_lam_vec)
% [v,v_dot] = make_arbitrary_steering_vector(d_lam_vec,N)
%
% Returns an array manifold for an arbitrary linear array with 
% inter-element spacing defined by the vector d_lam_vec.
%
% Inputs:
%
%   d_lam_vec   Inter-element spacing, in units of wavelengths, as a vector
%               of length N-1 for an N element array
%
% Outputs:
%
%   v           Function handle that accepts an angle (in radians) and
%               returns an N-element vector of complex phase shifts for
%               each element.  If multiple angles are supplied, the output
%               is a matrix of size numel(th) x N.
%   v_dot       Function handle to gradient vector
%               dv(psi) / dpsi
%
% Nicholas O'Donoughue
% 1 July 2019

% Ensure that the input is a row vector
d_lam_vec = reshape(d_lam_vec,1,numel(d_lam_vec));

% Accumulate the inter-element spacings, to compute the distance from each
% element to the reference element (first)
cum_spacing = [0 cumsum(d_lam_vec)];

v = @(psi) exp(1i*2*pi*cum_spacing.*sin(psi(:)));
v_dot = @(psi) (1i*2*pi*cum_spacing.*cos(psi(:))).*v(psi);