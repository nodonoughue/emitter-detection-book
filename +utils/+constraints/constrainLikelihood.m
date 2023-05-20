function ellC = constrainLikelihood(ell, a, b, tol)
%ellC = constrainLikelihood(ell, a, b, tol)
%
% Accepts a set of functions handles ell (likelihood), a (equality 
% constraints), b (inequality constraints), a tolerance to apply
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

%% Parse Inputs

% Check for an equality constraint
doEqConst = nargin >= 2 && ~isempty(a);

% Check for an inequality constraint
doIneqConst = nargin >= 3 && ~isempty(b);

% If there is an equality constraint, we need a defined tolerance
assert(~doEqConst || nargin >= 4 && ~isempty(tol), ...
       'No tolerance provided; equality constraints require a tolerance.');

% Initialize the constraints/tolerances, if not provided.
if ~doEqConst, a = []; tol = []; end
if ~doIneqConst, b = []; end

ellC = @(x) constrainLikelihoodInner(ell, x, a, b, tol, doEqConst, doIneqConst);

end
%% Construct Likelihood Constraint Function and Return a Handle
function result = constrainLikelihoodInner(ell, x, a, b, tol, doEqConst, doIneqConst)
    % Apply Constraints
    if doEqConst
        mask = checkEqConst(x, a, tol);
    else
        mask = 1;
    end

    if doIneqConst
        mask = mask & checkIneqConst(x, b);
    end

    if numel(mask) == 1
        % Something went wrong, just compute the likelihood
        result = ell(x);
        return;
    end
    
    % Initialize output; overwrite (for points that are valid)
    result = -Inf * ones(size(x));
    result(mask) = ell(x(:,mask));
end

%% Define Function to Check Equality Constraint
function eq_mask = checkEqConst(x, a, tol)
    % Check the Equality Constraint, and return a binary mask -- true for
    % all points that meet the constraint and false for all those that fail
    eq_mask = true(size(x,2),1);
    if numel(a) == 1, a = {a}; end

    for idxConst = 1:numel(a)
        % Grab the current constraint (function handle)
        const = a{idxConst};

        % Test only those points that are still valid (haven't violated any
        % previous constraints
        [eps, ~] = const(x(:,eq_mask));
        this_valid_mask = abs(eps) <= tol;
        
        % Overwrite the set of mask points that are currently valid with
        % the outcome of this constraint test
        eq_mask(eq_mask) = this_valid_mask;
    end
end

%% Define Function to Check Inequality Constraint
function ineq_mask = checkIneqConst(x, b)   
    % Check the Inequality Constraint, and return a binary mask -- true for
    % all points that meet the constraint and false for all those that fail
    ineq_mask = true(size(x,2),1);
    
    if numel(b) == 1, b = {b}; end

    % Check all constraints and mask out any points that violate them.
    for idxConst = 1:numel(b)
        % Grab the current constraint (function handle)
        const = b{idxConst};

        % Test only those points that are still valid (haven't violated any
        % previous constraints
        [eps, ~] = const(x(:,ineq_mask));
        this_valid_mask = eps > 0;
        
        % Overwrite the set of mask points that are currently valid with
        % the outcome of this constraint test
        ineq_mask(ineq_mask) = this_valid_mask;
    end
end