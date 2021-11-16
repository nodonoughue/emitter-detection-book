function [x,x_full] = lsSolnFixed(x_fdoa,v_fdoa,rho_dot,C,x_init,a,tol,epsilon,max_num_iterations,force_full_calc,plot_progress,ref_idx)
% [x,x_full] = lsSolnFixed(x_fdoa,v_fdoa,rho_dot,C,x_init,epsilon,numIterations,
%                                 force_full_calc,plot_progress,ref_idx)
%
% Computes the least square solution for FDOA processing.
%
% Utilized the utils.constraints package to accept equality constraints (a)
% with tolerance (tol).
%
% Inputs:
%   
%   x_fdoa              Sensor positions [m]
%   v_fdoa              Sensor velocities [m/s]
%   rho_dot             Range Rate-Difference Measurements [m/s]
%   C                   Measurement Error Covariance Matrix [(m/s)^2]
%   x_init              Initial estimate of source position [m]
%   a               Array of equality constraint function handles
%   tol             Tolerance for equality constraints
%   epsilon             Desired estimate resolution [m]
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   ref_idx             Scalar index of reference sensor, or nDim x nPair
%                       matrix of sensor pairings  
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 14 Nov 2021

% Parse inputs
if nargin < 12 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 11 || ~exist('plot_progress','var')
    plot_progress = false;
end

if nargin < 10 || ~exist('force_full_calc','var')
    force_full_calc = false;
end

if nargin < 8 || ~exist('max_num_iterations','var')
    max_num_iterations = [];
end

if nargin < 8 || ~exist('epsilon','var')
    epsilon = [];
end

if nargin < 7 || ~exist('tol','var')
    tol = [];
end

% Initialize measurement error and Jacobian function handles
y = @(x) rho_dot - fdoa.measurement(x_fdoa, v_fdoa, x, ref_idx);
J = @(x) fdoa.jacobian(x_fdoa, v_fdoa, x, ref_idx);

% Resample covariance matrix
n_sensor = size(x_fdoa, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Call the generic Least Square solver
[x,x_full] = utils.constraints.lsSolnFixed(y,J,C_tilde,x_init,a, tol,epsilon,max_num_iterations,force_full_calc,plot_progress);