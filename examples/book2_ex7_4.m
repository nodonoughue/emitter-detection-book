function figs = book2_ex7_4()
% figs = book2_ex7_4()
%
% Executes Example 7.4 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 11 January 2022

fprintf('Example 7.4...\n');

% Define sensor positions
x_aoa = [-5e3, 5e3;
         0, 0];
v_aoa = [0, 0;
         200, 200];
x_tgt = [ 3e3;
         25e3];

% Define sensor accuracy
num_dims = size(x_aoa,1);
num_sensors = size(x_aoa,2);
sigma_theta = 10;
sigma_psi = sigma_theta*pi/180;
C_df = sigma_psi^2*eye(num_sensors);
L = chol(C_df,'lower');

% Define pulse timing
pri = 1e-3;
T = 30; % observation period
num_pulses = floor(T/pri)+1;

% Define outputs
x_est = zeros(size(x_tgt,1),num_pulses);
cep = zeros(1,num_pulses);
x_init = [1; 1]*1e3;

% Step through pulses
for idx=1:num_pulses
    % Update positions
    this_x_aoa = x_aoa + v_aoa * (idx-1) * pri;

    % Update function handles
    z_fun = @(x) triang.measurement(this_x_aoa, x);
    H_fun = @(x) triang.jacobian(this_x_aoa, x)';

    % Generate noisy measurements
    psi = z_fun(x_tgt);
    zeta = psi + L*randn(num_dims,1);

    if idx==1
        % Initialization
        this_x = triang.lsSoln(this_x_aoa, zeta, C_df, x_init);
        this_P = triang.computeCRLB(this_x_aoa, this_x, C_df);
    else
        % EKF Update
        [this_x, this_P] = tracker.ekfUpdate(prev_x, prev_P, zeta, C_df, z_fun, H_fun);
    end

    % Store the results and update the variables
    x_est(:, idx) = this_x;
    cep(idx) = utils.computeCEP50(this_P);

    prev_x = this_x;
    prev_P = this_P;
end

fig1=figure;
plot(x_est(1,:),x_est(2,:),'-.','DisplayName','Estimated Position');
hold on;
plot(x_tgt(1),x_tgt(2),'^','DisplayName','Target');
grid on;
legend('Location','NorthEast');

% Draw Error Ellipse from single sample
crlb = triang.computeCRLB(this_x_aoa, x_tgt, C_df);
ell = utils.drawErrorEllipse(x_tgt, crlb, 101);
plot(ell(1,:), ell(2,:),'-.','DisplayName','Error Ellipse (single msmt.)');

ell_end = utils.drawErrorEllipse(x_tgt, this_P, 101);
plot(ell_end(1,:), ell_end(2,:),'-.','DisplayName','Error Ellipse (Final EKF Update)');

offset = max(max(ell,[],2)-min(ell,[],2));
xlim(x_tgt(1) + .6*offset*[-1 1]);
ylim(x_tgt(2) + .6*offset*[-1 1]);
utils.setPlotStyle(gca,{'widescreen'});

%% Compute Errors
err = sqrt(sum(abs(x_est - x_tgt).^2,1));
time_vec = pri*(1:num_pulses);

fig2=figure;
semilogy(time_vec,err,'DisplayName','Measured');
hold on;
plot(time_vec,cep,'DisplayName','Predicted (CEP_{50})');
xlabel('Time [s]');
ylabel('Error [m]');
legend('Location','NorthEast');
grid on;
utils.setPlotStyle(gca,{'widescreen'});

%% Find Time to Achieved Desired Error
desired_cep = 100;
first_good_sample = find(cep < desired_cep,1,'first');
if isempty(first_good_sample)
    fprintf('More than %.2f s required to achieve %.2f m CEP50.\n', max(time_vec), desired_cep);
else
    fprintf('%.2f s required to achieve %.2f m CEP50.\n', time_vec(first_good_sample), desired_cep);
end


figs = [fig1 fig2];