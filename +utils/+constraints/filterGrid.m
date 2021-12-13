function [x_filtered, valid_mask] = filterGrid(x_grid, a, b, tol)
%FILTERGRID Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2 || isempty(a)
    doEqConst = false;
    a = @(x) 0;
else
    doEqConst = true;
end

if nargin < 3 || isempty(b)
    doIneqConst = false;
    b = @(x) 0;
else
    doIneqConst = true;
end

if (nargin < 4 || isempty(tol)) && doEqConst
    % No tolerance provided; compute one based on average grid spacing
    dist = sqrt(sum(abs(x_grid(:,2:end) - x_grid(:,1:end-1)).^2,1));
    min_dist = min(dist);
   
    tol = sqrt(2) * min_dist;
end

% Initialize output mask
n_points = size(x_grid,2);
n_dim = size(x_grid,1);
valid_mask = true(1,n_points);

% Equality Constraint Testing
if doEqConst
    for idx = 1:n_points
        for idxConst = 1:numel(a)
            const = a(idxConst);
            [eps, ~] = const(x_grid(:,idx));
            if abs(eps) > tol
                valid_mask(idx) = false; 
                break;
            end
        end
    end
end

if doIneqConst
   for idx = 1:n_points
       if ~valid_mask(idx), continue, end
       for idxConst = 1:numel(b)
           const = b(idxConst);
           [eps, ~] = const(x_grid(:,idx));
           if eps > 0
               valid_mask(idx) = false;
               break;
           end
       end
   end
end


end

