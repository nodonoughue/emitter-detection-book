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

if doEqConst && ~doIneqConst
    ellC = @(x) constrainedLikelihoodEq(x, ell, a, tol);
elseif doEqConst
    ellC = @(x) constrainedLikelihoodEqIneq(x, ell, a, tol, b);
else
    ellC = @(x) constrainedLikelihoodIneq(x, ell, b);
end

return


function ell_out = constrainedLikelihoodEq(x, ell, a, tol)
    % Start with unconstrained likelihood
    ell_out = ell(x);
    
    % Check all constraints and set likelihood to -Inf for any point
    % that violate them.
    for idxConst = 1:numel(a)
        const = a(idxConst);
        [eps, ~] = const(x);
        valid_mask = abs(eps) <= tol;
        ell_out(~valid_mask) = -Inf;
    end
        
end

function ell_out = constrainedLikelihoodEqIneq(x, ell, a, tol, b)
    % Start with unconstrained likelihood
    ell_out = ell(x);
    
    % Check all eq constraints and set likelihood to -Inf for any point
    % that violate them.
    for idxConst = 1:numel(a)
        const = a(idxConst);
        [eps, ~] = const(x);
        valid_mask = abs(eps) <= tol;
        ell_out(~valid_mask) = -Inf;
    end

    % Check all ineq constraints and set likelihood to -Inf for any point
    % that violate them.
    for idxConst = 1:numel(b)
        const = b(idxConst);
        [eps, ~] = const(x);
        valid_mask = eps <= 0;
        ell_out(~valid_mask) = -Inf;
    end

end

function ell_out = constrainedLikelihoodIneq(x, ell, b)
    % Start with unconstrained likelihood
    ell_out = ell(x);
    
    % Check all ineq constraints and set likelihood to -Inf for any point
    % that violate them.
    for idxConst = 1:numel(b)
        const = b(idxConst);
        [eps, ~] = const(x);
        valid_mask = eps <= 0;
        ell_out(~valid_mask) = -Inf;
    end
end

end