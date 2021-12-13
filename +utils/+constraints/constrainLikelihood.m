function ellC = constrainLikelihood(ell, a, b, tol)
%ellC = constrainLikelihood(ell, a, b, tol)
%
% Accepts a set of functions handles ell (likelihood), a (equality 
% constraints), b (inequality constraints), and a tolerance to apply
% to the equality constraints.
%
% Returns a contrained likelihood function handle that will accept an
% nDim x N set of position vectors, and return a vector of N constrained
% likelihood values.  If either the abs(a(x))<= tol, or b(x)<=0 constraints
% are violated, then the likelihood is -Inf.
%
% INPUTS:
%   ell     Likelihood function handle
%   a       Equality constraint function handle
%   b       Inequality constraint function handle
%   tol     Tolerance for equality constraint
%
% Nicholas O'Donoughue
% 5 Sept 2021


if nargin < 2 || isempty(a)
    doEqConst = false;
    a = @(x) zeros(1,size(x,2));
else
    doEqConst = true;
end

if nargin < 3 || isempty(b)
    doIneqConst = false;
    b = @(x) zeros(1,size(x,2));
else
    doIneqConst = true;
end

if (nargin < 4 || isempty(tol)) && doEqConst
    error('No tolerance provided; equality constraints require a tolerance.');
end

%% Determine which points are valid
valid_mask = true(size(x));

if doEqConst
    % Check all constraints and mask out any points that violate them.
    for idxConst = 1:numel(a)
        % Grab the current constraint (function handle)
        const = a{idxConst};

        % Test only those points that are still valid (haven't violated any
        % previous constraints
        [eps, ~] = const(x(valid_mask));
        this_valid_mask = abs(eps) <= tol;
        
        % Overwrite the set of mask points that are currently valid with
        % the outcome of this constraint test
        valid_mask(valid_mask) = this_valid_mask;
    end
end

if doIneqConst
    % Check all constraints and mask out any points that violate them.
    for idxConst = 1:numel(b)
        % Grab the current constraint (function handle)
        const = b{idxConst};

        % Test only those points that are still valid (haven't violated any
        % previous constraints
        [eps, ~] = const(x(valid_mask));
        this_valid_mask = eps > 0;
        
        % Overwrite the set of mask points that are currently valid with
        % the outcome of this constraint test
        valid_mask(valid_mask) = this_valid_mask;
    end
end

%% Evaluate the Likelihood
ellC = -inf * ones(size(x));
ellC(valid_mask) = ell(x(valid_mask));        
