function [P,psi_vec] = music(x,v,D,max_psi,N_pts)
% [P,psi_vec] = music(x,v,D,max_psi,N_pts)
%
% Generates a MUSIC-based image for N_pts equally spaced angular 
% coordinates from -psi_max to psi_max (in radians), given the input 
% data x, array steering vector v, and optional number of signals D.
%
% If left blank, or set to zero, the value D will be estimated using a
% simple algorithm that counts the number of eigenvalues greater than twice
% the minimum eigenvalue.  This will break down in low SNR scenarios.
%
% INPUTS:
%
%   x       N x M data vector
%   v       Steering vector function that returns
%           N point steering vector for each input (in radians).
%   D       Number of signals [optional, set to zero to automatically
%           estimate D based on a threshold eigenvalue twice the minimum
%           eigenvalue]
%   max_psi Maximum steering angle (radians)
%   N_pts   Number of steering angles to compute
%
% OUTPUTS:
%
%   P       Power image (1 x N_pts) in linear units
%   psi_vec Vector of scan angles computed (in radians)
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 5 || isempty(N_pts)
    N_pts = 101;
end

if nargin < 4 || isempty(max_psi)
    max_psi = pi/2;
end

if nargin < 3 || isempty(D)
    D = 0;
end

% Compute the sample covariance matrix
[N,M] = size(x);
C = zeros(N,N);
for idx_m = 1:M
    C = C + x(:,idx_m)*x(:,idx_m)';
end
C = C/M;

% Perform Eigendecomposition of C
[U,Lam] = eig(C);
lam = diag(Lam);

% Sort the eigenvalues
[~,idx_sort] = sort(abs(lam),'descend');
U_sort = U(:,idx_sort);
lam_sort = lam(idx_sort);

% Isolate Noise Subspace
if D ~=0
    Un = U_sort(:,D+1:end);
else
    % We need to estimate D first
    
    % Assume that the noise power is given by the smallest eigenvalue
    noise = min(lam);
    
    % Set a threshold of 2x the noise level; and find the eigenvalue that
    % first cuts above it
    D = find(lam_sort >= 2*noise,1,'last');
    
    Un = U_sort(:,D+1:end);
end

% MUSIC Projection
proj = Un*Un';

% Generate steering vectors
psi_vec = linspace(-1,1,N_pts)*max_psi;
P = zeros(1,N_pts);
for idx_pt = 1:N_pts
    vv = v(psi_vec(idx_pt))/sqrt(N);
    Q = vv'*proj*vv;
    P(idx_pt) = 1./abs(Q);
end
