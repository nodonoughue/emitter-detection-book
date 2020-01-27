function psi = directional_df(s,psi_samples,g,psi_res,psi_min,psi_max)
% psi = directional_df(s,psi_samples,g,psi_res,psi_min,psi_max)
%
% Computes an estimate of the angle of arrival psi (in radians) for a set
% of amplitude measurements s, taken at various steering angles psi_samples
%
% Inputs:
%
%   s           Set of M measurements taken at each of N steering angles,
%               contains M rows, and N columsn.
%   psi_samples Steering angles at which measurements were taken [radians]
%   g           Function handle to gain equation, g(psi)
%   psi_res     Desired resolution for output estimate [defualt = .1]
%   psi_min     Lower bound on valid region for psi [default = -pi]
%   psi_max     Upper bound on valid region for psi [default = pi]
%
% Outputs:
%
%   psi         Estimated angle of arrival [radians]
%
% Nicholas O'Donoughue
% 1 July 2019

% Handle incomplete inputs
if nargin < 6 || isempty(psi_max)
    psi_max = pi;
end
if nargin < 5 || isempty(psi_min)
    psi_min = -pi;
end
if nargin < 4 || isempty(psi_res)
    psi_res = .1;
end

% Determine how many samples exist
[N,M] = size(s);

% Initialize the outer loop to .1 radians
this_psi_res = .1;

% Iteratively search for the optimal index, until the sample spacing is
% less than the desired resolution
while this_psi_res >= psi_res/10
    % Define the vector of possible AoA samples for the current resolution
    psi_vec = psi_min:this_psi_res:psi_max;

    % Find the difference between each possible AoA (psi_vec) and the
    % measurement steering angles
    psi_diff = psi_samples(:) - psi_vec;
    
    % Compute the expected gain pattern for each candidate AoA value
    g_vec = arrayfun(@(x) g(x),psi_diff);

    % Find the candidate AoA value that minimizes the MSE between the
    % expected and received gain patterns.
    mse = squeeze(sum(sum(abs(reshape(s,[N,1,M])-g_vec).^2,1),3));
    [~,idx] = min(mse);
    psi=psi_vec(idx);
    
    % Set up the bounds and resolution for the next iteration
    psi_min = psi_vec(max(1,idx-4));
    psi_max = psi_vec(min(numel(psi_vec),idx+4));
    this_psi_res = this_psi_res/10;
end