function x_valid = snapToEqConstraint(x, a, tol)

for idx=1:numel(a)
    this_a = a(idx);
    
    [epsilon, scale] = this_a(x);
    invalid_mask = abs(epsilon) > tol;
    
    x(invalid_mask) = x(invalid_mask)*scale(invalid_mask);
end

x_valid = x;