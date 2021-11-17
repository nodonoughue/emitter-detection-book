function b = boundedAlt(h_min, h_max, type)

% Upper bound
if ~isempty(h_max) && isfinite(h_max)
    [a_upper, ~] = utils.constraints.fixedAlt(h_max, type);

    b = {a_upper};
else
    b = [];
end

% Lower bound
if ~isempty(h_min) && isfinite(h_min)
    [a_lower, ~] = utils.constraints.fixedAlt(h_min, type);

    % Invert the epsilon sign on the lower bound (so that <= 0 is what we want)
    a_lower_inv = @(x) invertEpsilon(a_lower, x);

    % Return the full set of inequality bounds
    b = cat(1,b, {a_lower_inv});
end

function [epsilon_inverse, x_valid] = invertEpsilon(a, x)

    [epsilon, x_valid] = a(x);
    epsilon_inverse = -epsilon;

end

end