function figs = book2_ex7_2()
% figs = book2_ex7_2()
%
% Executes Example 7.2 from Practical Geolocation for Electronic Warfare
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

fprintf('Example 7.2...\n');

% Define sensor positions
x_tdoa = [-10e3, 0, 10e3;
         0, 10e3, 0];
x_tgt = [  5e3;
         -15e3];

% Define sensor accuracy
n_sensors = size(x_tdoa,2);
sigma_t = 1e-6;
C_toa = sigma_t^2*eye(n_sensors);
ref_idx = []; % use default reference sensor

% Define pulse timing
pri = 1e-3;

% Compute CRLB
crlb_single_sample = tdoa.computeCRLB(x_tdoa, x_tgt, C_toa, ref_idx);
cep_single_sample = utils.computeCEP50(crlb_single_sample);
fprintf('CEP50 for a single sample: %.2f km\n',cep_single_sample/1e3)

% Iterate over observation interval
time_vec = [.1:.1:100, 125:25:1000]; % seconds
sigma_t_vec = [.1e-6, 1e-6, 10e-6];
cep_vec = zeros(numel(time_vec), numel(sigma_t_vec));
for idx_t = 1:numel(time_vec)
    this_time = time_vec(idx_t);
    this_num_samples = 1 + floor(this_time/pri);
    for idx_s = 1:numel(sigma_t_vec)
        this_sigma_t = sigma_t_vec(idx_s);
        this_C_toa = this_sigma_t^2*eye(n_sensors);
        this_C = this_C_toa/this_num_samples;

        this_crlb = tdoa.computeCRLB(x_tdoa, x_tgt, this_C, ref_idx);
        cep_vec(idx_t, idx_s) = utils.computeCEP50(this_crlb);
    end
end

fig1=figure;
for idx_s = 1:numel(sigma_t_vec)
    loglog(time_vec, cep_vec(:,idx_s),'DisplayName',...
             sprintf('\\sigma_t=%.1f \\mu s',sigma_t_vec(idx_s)*1e6));
    hold on;
end
plot(time_vec, 10*ones(size(time_vec)),'k-.','DisplayName','10 m');
xlabel('Time [s]');
ylabel('$CEP_{50}$ [m]');
ylim([1, 1000]);
legend('Location','NorthEast');
grid on;
utils.setPlotStyle(gca,{'widescreen'});

% Determine when CEP50 crosses below 10 m
desired_cep = 10;
first_good_sample = find(cep_vec <= desired_cep, 1, 'first');
if isempty(first_good_sample)
    fprintf('More than %.2f s required to achieve %.2f m CEP50.\n', max(time_vec), desired_cep);
else
    fprintf('%.2f s required to achieve %.2f m CEP50.\n', time_vec(first_good_sample), desired_cep);
end

T = time_vec(first_good_sample);
num_pulses = floor(T/pri)+1;

% Compute CRLB
crlb_sample_mean = tdoa.computeCRLB(x_tdoa, x_tgt, C_toa/num_pulses, ref_idx);
cep_sample_mean = utils.computeCEP50(crlb_sample_mean);
fprintf('CEP50 for the sample mean: %.2f km\n',cep_sample_mean/1e3);


%% Demonstrate geolocation
z = tdoa.measurement(x_tdoa, x_tgt, ref_idx);
C_roa = utils.constants.c^2*C_toa;
C_rdoa = utils.resampleCovMtx(C_roa,[]);
L = chol(C_rdoa,'lower');
zeta = z + L*randn(size(C_rdoa,1),num_pulses); % noisy measurement

% Sample Mean
zeta_mn = cumsum(zeta,2)./(1:num_pulses);

% Geolocation Result
x_ls = zeros(numel(x_tgt),num_pulses);
x_ls_mn = zeros(numel(x_tgt), num_pulses);

x_init = [1; 1]*1e3;

for idx=1:num_pulses
    [x_ls(:,idx), ~] = tdoa.lsSoln(x_tdoa, zeta(:,idx), C_roa, x_init);
    [x_ls_mn(:,idx), ~] = tdoa.lsSoln(x_tdoa, zeta_mn(:,idx), C_roa/idx, x_init);
end

fig2=figure;
% plot(x_tdoa(1,:), x_tdoa(2,:),'o','DisplayName','TDOA Sensors');
% hold on;
plot(x_tgt(1,:), x_tgt(2,:), '^','DisplayName','Target');
hold on;
%plot(x_ls_mn(1,:),x_ls_mn(2,:),'--','DisplayName','LS Soln (single sample)');
plot(x_ls_mn(1,:),x_ls_mn(2,:),'-.','DisplayName','LS Soln (sample mean)');
grid on;
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen','equal'});

% Overlay error ellipse
ell = utils.drawErrorEllipse(x_tgt,crlb_single_sample,101);
ell_full = utils.drawErrorEllipse(x_tgt,crlb_sample_mean,101);

plot(ell(1,:),ell(2,:),'DisplayName','Error Ellipse (single sample)');
plot(ell_full(1,:), ell_full(2,:),'DisplayName','Error Ellipse (sample mean)');

% Plot error as a function of time
err = sqrt(sum(abs(x_ls-x_tgt).^2,1));
err_mn = sqrt(sum(abs(x_ls_mn-x_tgt).^2,1));

fig3=figure;
time_vec = pri*(1:num_pulses);
semilogy(time_vec,err,'DisplayName','Error (single sample)');
hold on;
semilogy(time_vec,err_mn,'DisplayName','Error (sample mean)');
plot(time_vec,cep_single_sample*ones(1,num_pulses),'DisplayName','CRLB (single sample)');
plot(time_vec,cep_single_sample./sqrt(1:num_pulses),'DisplayName','CRLB (sample mean)');
xlabel('Time [s]');
ylabel('Error [m]');
ylim([10, 10e3]);
legend('Location','NorthEast');
grid on;
utils.setPlotStyle(gca,{'widescreen'});

figs = [fig1, fig2, fig3];
return