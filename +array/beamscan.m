function [P,psi_vec] = beamscan(x,v,psi_max,N_pts)
% [P,psi_vec] = beamscan(x,v,psi_max,N_pts)
%
% Generates a beamscan image for N_pts equally spaced angular coordinates
% from -psi_max to psi_max (in radians), given the input data x, and array
% steering vector v
%
% INPUTS:
%
%   x       N x M data vector
%   v       Steering vector function that returns
%           N point steering vector for each input (in radians).
%   psi_max Maximum steering angle (radians)
%   N_pts   Number of steering angles to compute
%
% OUTPUTS:
%
%   P       Power image (1 x N_pts) in linear units
%   psi_vec Vector of scan angles computed (in radians)
%
% Nicholas O'Donoughue
% 1 July 2019

% Generate scan vector
if nargin < 4 || isempty(N_pts)
    N_pts = 101;
end

if nargin < 3 || isempty(psi_max)
    psi_max = pi/2;
end

psi_vec = linspace(-1,1,N_pts)*psi_max;

% Parse inputs
[N,M] = size(x);

% Generate steering vectors
V = v(psi_vec)/sqrt(N); % N x N_pts

% Steer each of the M data samples
% - take the magnitude squared
% - compute the mean across M snapshots
P = sum(abs(x'*V).^2,1)/M; % 1 x N_pts
