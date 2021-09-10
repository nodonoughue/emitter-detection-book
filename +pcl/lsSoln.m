function [x,x_full] = lsSoln(x_tx,x_rx,rho,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,ref_idx)
% [x,x_full] = lsSoln(x_tdoa,rho,C,x_init,epsilon,max_num_iterations,...
%                               force_full_calc, plot_progress,ref_idx)
%
% Computes the least square solution for PCL processing -- range
% measurements only for now.
%
% TODO: Incorporate range rate measurements
% TODO: Incorporate angle of arrival measurements
%
% Inputs:
%   
%   x_tx                Transmitter positions [m]
%   x_rx                Receiver positions [m]
%   rho                 Bistatic Range Measurements [m]
%   C                   Measurement Error Covariance Matrix [m^2]
%   x_init              Initial source position estimate [m]
%   epsilon             Desired estimate resolution [m]
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   ref_idx             Matrix of tx/rx pairing indices (in the order
%                       they're used in C).  If ignored, then all pairwise
%                       measurements are used (nTx x nRx)
% 
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 10 September 2021

% Parse inputs
if nargin < 10 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 9 || ~exist('plot_progress','var')
    plot_progress = false;
end

if nargin < 8 || ~exist('force_full_calc','var')
    force_full_calc = false;
end

if nargin < 7 || ~exist('max_num_iterations','var')
    max_num_iterations = [];
end

if nargin < 6 || ~exist('epsilon','var')
    epsilon = [];
end

% Initialize measurement error and Jacobian function handles
y = @(x) rho - pcl.measurement(x_tx, x_rx, x, ref_idx);
J = @(x) pcl.jacobian(x_tx, x_rx, x, [], [], [], ref_idx);

% Call the generic Least Square solver
[x,x_full] = utils.lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);