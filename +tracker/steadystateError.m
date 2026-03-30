function [P,Pe] = steadystateError(F,H,Q,R,max_num_iterations,epsilon)
% [P,Pe] = steadystateError(F,H,Q,R,max_num_iterations)
%
% Iteratively compute the steady-state prediction error for a Kalman Filter
% based on the steady-state Ricatti equation for prediction error:
%       P = Q + F( P^{-1} + H'R^{-1}H )^{-1} F'
%
% where
%   F = state transition matrix
%   H = measurement matrix
%   Q = process noise (variance of state from one update to the next)
%   R = measurement noise (variance of each set of measurements y)
%
% This is done iteratively, until the difference between successive
% estimates of F is sufficiently small, or until the max_num_iterations is
% reached.
%
% INPUTS:
%   F       N x N state transition matrix
%   H       M x N measurement matrix
%   Q       N x N covariance of process noise
%   R       M x M covariance of measurement noise
%   max_num_iterations (optional) stopping condition for iteration
%   epsilon (optional) stopping condition for iteration
%
% OUTPUTS:
%   P       N x N steady-state Prediction error covariance
%   Pe      N x N steady-state Estimation error covariance
%
% Nicholas O'Donoughue
% 11 Nov 2021

%% Parse Inputs
if nargin < 6 || ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 1e-6 * norm(F,'fro');
end

if nargin < 5 || ~exist('max_num_iterations','var') || isempty(max_num_iterations)
    max_num_iterations = 100;
end

isOctave = exist('OCTAVE_VERSION','builtin')~=0;

%% Check Input Dimensions
%  Allow, F, H, Q, and R to have arbitrary size.  The first two dimensions
%  are used for the tracker.  Any extra non-singleton dimensions must
%  match.

f_dims = size(F);
h_dims = size(H);
q_dims = size(Q);
r_dims = size(R);

if numel(f_dims) > 2 || numel(h_dims) > 2 || numel(q_dims) > 2 || numel(r_dims) > 2
    % At least one input has >2 dimensions
    if numel(f_dims) > 2
        f_case_dims = f_dims(3:end);
    else
        f_case_dims = 1;
    end
    
    if numel(h_dims) > 2
        h_case_dims = h_dims(3:end);
    else
        h_case_dims = 1;
    end
    
    if numel(q_dims) > 2
        q_case_dims = q_dims(3:end);
    else
        q_case_dims = 1;
    end
    
    if numel(r_dims) > 2
        r_case_dims = r_dims(3:end);
    else
        r_case_dims = 1;
    end
    
    % Pad the case dimensions
    n_case_dims = max([numel(f_case_dims),numel(h_case_dims),numel(q_case_dims),numel(r_case_dims)]);
    f_case_dims = [f_case_dims, ones(1,n_case_dims-numel(f_case_dims))];
    h_case_dims = [h_case_dims, ones(1,n_case_dims-numel(h_case_dims))];
    q_case_dims = [q_case_dims, ones(1,n_case_dims-numel(q_case_dims))];
    r_case_dims = [r_case_dims, ones(1,n_case_dims-numel(r_case_dims))];
    
    % Find the cumulative case dims
    case_dims = max([f_case_dims; h_case_dims; q_case_dims; r_case_dims],[],1);
    
    % Check compatibility
    assert(all(f_case_dims == 1 | f_case_dims==case_dims) && ...
           all(h_case_dims == 1 | h_case_dims==case_dims) && ...
           all(q_case_dims == 1 | q_case_dims==case_dims) && ...
           all(r_case_dims == 1 | r_case_dims==case_dims), ...
           'All non-singleton dimensions (after the first two) of Kalman Filter matrices must match.');
       
    % Pad the inputs
    F = repmat(F,[1,1,case_dims./f_case_dims]);
    H = repmat(H,[1,1,case_dims./h_case_dims]);
    Q = repmat(Q,[1,1,case_dims./q_case_dims]);
    R = repmat(R,[1,1,case_dims./r_case_dims]);
    
    % Use arrayfun to execute
    out_dims = [f_dims(1:2), case_dims];
    P = zeros(out_dims);
    Pe = zeros(out_dims);
    for idx_case = 1:prod(case_dims)
        % Use recursion to loop over cases
        [thisP, thisPe] = tracker.steadystateError(F(:,:,idx_case),H(:,:,idx_case),Q(:,:,idx_case),R(:,:,idx_case),max_num_iterations,epsilon);
        P(:,:,idx_case) = thisP;
        Pe(:,:,idx_case) = thisPe;
    end
    return;
end


%% Initialize P
P = eye(size(Q));
dist = inf;
iter = 0;

%% Pre-compute R inverse
HRH = (H' / R) * H;

%% Loop until convergence
while iter < max_num_iterations && abs(dist) > epsilon
    % Update estimate
    P_prev = P;
    
    % Pre-invert P_prev, using knowledge that it's a symmetric matrix
    P = Q + (F / (inv(P_prev) + HRH)) * F';
    
    % Ensure P stays symmetric
    P = .5 * (P + P.');
    
    % Increment the convergence counters
    dist = norm(P-P_prev,'fro');
    iter = iter + 1;
end

%% Comput Estimation Error
Pe = inv(inv(P) + HRH);