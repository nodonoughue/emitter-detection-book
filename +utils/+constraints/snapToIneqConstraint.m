function x_valid = snapToIneqConstraint(x, b)

for idx=1:numel(b)
    this_b = b(idx);
    
    [epsilon, scale] = this_b(x);
    invalid_mask = epsilon > 0;
    
    x(invalid_mask) = x(invalid_mask)*scale(invalid_mask);
end

x_valid = x;