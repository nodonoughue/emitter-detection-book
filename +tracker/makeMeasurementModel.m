function [z_fun, h_fun] = makeMeasurementModel(x_aoa, x_tdoa, x_fdoa, v_fdoa, tdoa_ref_idx, fdoa_ref_idx, state_space)

% Sample the position/velocity components of the target state
x_pos = @(x) x(state_space.pos_idx,:);
if state_space.has_vel
    x_vel = @(x) x(state_space.vel_idx,:);
else
    x_vel = @(x) [];
end

if isempty(v_fdoa)
    v_fdoa = zeros(numel(state_space.vel_idx),1);
end

%% Non-Linear Measurement Function
z_fun = @(x) hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa - x_vel(x), x_pos(x), tdoa_ref_idx, fdoa_ref_idx);
    % num_msmt x size(x,2)

%% Measurement Function Generator
function h = buildMsmtFunction(x)
    [J, Jv] = hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_pos(x), tdoa_ref_idx, fdoa_ref_idx, x_vel(x));
    % num_dims x num_msmt 
    % The first is the gradient w.r.t. position, the second is w.r.t.
    % velocity

    [~, num_msmt, num_src] = size(J);

    h = zeros(num_msmt, state_space.num_states, num_src);
    h(:, state_space.pos_idx,:) = J.';

    if ~isempty(Jv)
        h(:, state_space.vel_idx,:) = Jv.';
    end
end

% Return a function handle to the measurement function generator
h_fun = @buildMsmtFunction;

end   