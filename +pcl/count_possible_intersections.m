function n_comb = count_possible_intersections(n_targ, n_pair, n_inter)
% Computes the number of possible ellipsoid intersections to test for
% validity, given n_target true targets in the scene, n_pair sets of
% transmitter/receiver pairs (each generates an ellipsoid), and n_inter 
% intersections required for position estimation (n>=3 for SX method, >=4
% for SI method).
%
% Each target is assumed to generate a measurement at each n_pair
% transmitter/receiver pair.
%
% From Section 8.2 of Malanowski, "Signal Processing for Passive Bistatic
% Radar," Artech House, 2019.
%
% INPUTS:
%   n_targ          Number of true targets in the scene
%   n_pair          Number of transmitter/receiver pairs processed
%   n_inter         Number of ellipsoid intersections needed for a target
%                   solution (default = 3)
%
% OUTPUTS:
%   n_comb          Number of possible interactions that need to be
%                   tested/resolved
%
% Nicholas O'Donoughue
% 6 April 2020

if nargin < 3 || isempty(n_inter)
    n_inter = 3;
end

n_comb = n_targ.^n_inter .* arrayfun(@(x) nchoosek(x, n_inter),n_pair);