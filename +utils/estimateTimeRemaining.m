function t_rem = estimateTimeRemaining(t_elapsed, idx, N)
% Estimate time remaining, given elapsed time and progress.
%
%

if idx==0
    % Not enough info; impossible to estimate.
    t_rem = nan;
    return;
end

t_total = t_elapsed * N / idx;
t_rem = t_total - t_elapsed;