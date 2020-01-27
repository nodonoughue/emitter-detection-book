function [v,v_dot] = make_steering_vector(d_lam,N)
% [v,v_dot] = make_steering_vector(d_lam,N)
%
% Returns an array manifold for a uniform linear array with N elements
% and inter-element spacing d_lam.
%
% Inputs:
%
%   d_lam       Inter-element spacing, in units of wavelengths
%   N           Number of elements in array
%
% Outputs:
%
%   v           Function handle that accepts an angle psi (in radians) and
%               returns an N-element vector of complex phase shifts for
%               each element.  If multiple angles are supplied, the output
%               is a matrix of size N x numel(psi).
%   v_dot       Function handle that computes the gradient of v(psi) 
%               with respect to psi.  Returns a matrix of size N x
%               numel(psi)
%
% Nicholas O'Donoughue
% 1 July 2019

v = @(psi) exp(1i*2*pi*d_lam*(0:N-1).*sin(psi(:))).';

v_dot = @(psi) (1i*2*pi*d_lam*(0:N-1).*cos(psi(:))).'.*v(psi);