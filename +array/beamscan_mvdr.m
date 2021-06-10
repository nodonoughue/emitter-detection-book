function [P,psi_vec] = beamscan_mvdr(x,v,psi_max,N_pts)
% [P,psi_vec] = beamscan_mvdr(x,v,psi_max,N_pts)
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

% Compute the sample covariance matrix
[N,M] = size(x);
C = zeros(N,N);
for idx_m = 1:M
    C = C + x(:,idx_m)*x(:,idx_m)';
end
C = C/M;

% Make sure the covariance matrix is invertible
C = utils.ensureInvertible(C);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C);
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end

% Steer each of the M data samples
P = zeros(1,N_pts);
for idx_psi = 1:N_pts
    vv = v(psi_vec(idx_psi))/sqrt(N); % N x 1
    
    if do_decomp
        P(idx_psi) = 1./abs(vv'*(C_d\vv));
    else
        P(idx_psi) = 1./abs(vv'*C_inv*vv);
    end
end