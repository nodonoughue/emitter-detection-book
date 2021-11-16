function [x,x_full] = lsSolnFixed(y,J,C,x_init,a,tol,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = lsSolnFixed(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the least square solution for TDOA processing.
%
% Utilized the utils.constraints package to accept equality constraints (a)
% with tolerance (tol).
%
% Inputs:
%   
%   y               Measurement vector function handle (accepts n_dim 
%                   vector of source position estimate, responds with error 
%                   between received and modeled data vector)
%   J               Jacobian matrix function handle (accepts n_dim vector
%                   of source position estimate, and responds with n_dim x
%                   n_sensor Jacobian matrix
%   C               Measurement error covariance matrix
%   x_init          Initial estimate of source position
%   a               Array of equality constraint function handles
%   tol             Tolerance for equality constraints
%   epsilon         Desired position error tolerance (stopping condition)
%   max_num_iterations  Maximum number of LS iterations to perform
%   force_full_calc Forces all max_num_iterations to be calculated
%   plot_progress   Binary flag indicating whether to plot error/pos est
%                   over time
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 14 Nov 2021

% Parse inputs
n_dims = numel(x_init);

if nargin < 10 || isempty(plot_progress)
    % By default, do not plot
    plot_progress = false;
end

if nargin < 9 || isempty(force_full_calc)
    % If not specified, don't force the solver to run until it hits
    % max_num_iterations
    force_full_calc = false;
end

if nargin < 8 || isempty(max_num_iterations)
    % Number of iterations not specified, default is 10,000
    max_num_iterations = 10000;
end

if nargin < 7 || isempty(epsilon)
    % Maximum error not specified; default to .1 (distance units)
    epsilon = 1e-6;
end

% Initialize loop
iter = 1;
error = Inf;
x_full = zeros(n_dims,max_num_iterations);
x_prev = x_init;
x_full(:,1) = x_prev;
    
% For now, we assume the AOA is independent of TDOA/FDOA
C = utils.ensureInvertible(C);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end

% Initialize Plotting
if plot_progress
    figure;
    xlabel('Iteration Number');
    ylabel('Change in Position Estimate');
    hold on;
    set(gca,'yscale','log');
end

% Divergence Detection
num_expanding_iters = 0;
max_num_expanding_iters = 10;
prev_error = Inf;

% Loop until either the desired tolerance is achieved or the maximum
% number of iterations have occurred
while iter < max_num_iterations && (force_full_calc || error >= epsilon)
    iter = iter+1;

    % Evaluate Residual and Jacobian Matrix
    y_i = y(x_prev);
    J_i = J(x_prev);
            
    % Compute delta_x^(i), according to 13.18
    if do_decomp
        delta_x = (J_i/C_d*J_i')\(J_i/C_d)*y_i;
    else
        delta_x = (J_i*C_inv*J_i')\(J_i*C_inv)*y_i;
    end
    
    % Update predicted location
    x_unconst = x_prev + delta_x;

    % Apply Equality Constraints
    x_const = utils.constraints.snapToEqConstraint(x_unconst, a, tol);
    
    % Update variables
    x_full(:,iter) = x_const;
    x_prev = x_full(:,iter);
    error = norm(delta_x);

    if plot_progress
        plot(iter,error,'.');
    end
    
    % Check for divergence
    if error <= prev_error
        num_expanding_iters = 0;
    else
        num_expanding_iters = num_expanding_iters + 1;
        if num_expanding_iters >= max_num_expanding_iters
            % Divergence detected
            x_full(:,iter+1:end) = NaN;
            break;
        end
    end
    prev_error = error;
end

% Bookkeeping
if ~force_full_calc

    x_full = x_full(:,1:iter);
end
x = x_full(:,iter);
