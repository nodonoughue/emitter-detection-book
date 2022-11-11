% function figs = book2_ex7_1()
% figs=book2_ex7_1()
%
% Executes Example 7.1 from Practical Geolocation for Electronic Warfare
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

fprintf('Example 7.1...\n');
utils.resetPlotSettings; % reset the plot settings; they don't do well with subfigures

% Define sensor positions
x_aoa = [-2e3, 3e3;
         0, 0];
x_tgt = [2e3;
         4e3];

% Define sensor accuracy
n_dim = size(x_aoa,1);
n_aoa = size(x_aoa,2);
sigma_theta = 10;
sigma_psi = sigma_theta * pi/180;
C = sigma_psi^2 * eye(n_aoa);

% Compute CRLB
crlb = triang.computeCRLB(x_aoa, x_tgt, C);
cep = utils.computeCEP50(crlb)/1e3;

%% Initialize Variables
total_num_samples = 1000;

psi = triang.measurement(x_aoa,x_tgt);
noise = sqrt(C)*randn(n_aoa,total_num_samples);
zeta = psi + noise;

zeta_avg = cumsum(zeta,2)./(1:total_num_samples);

crlb_ellipses = arrayfun(@(k) utils.drawErrorEllipse(x_tgt, crlb/k,101),1:total_num_samples,'UniformOutput',false);
    % 1 x total_num_pulses cell array, eacho f which contains an error
    % ellipse

%% Solve Each One
x_prior = [0;1e3];
x_est_full = zeros(n_dim,total_num_samples);

for k=1:total_num_samples
    x_est = triang.lsSoln(x_aoa,zeta_avg(:,k),C,x_prior,[],[],false,false);
    x_est_full(:,k) = x_est;
%     x_prior = x_est;
end

%% Establish figure
fig=figure();
subplot(212);
xlabel('x [km]');
ylabel('y [km]');

% Plot the measurements
subplot(211);
hold on;
xlabel('Sample [K]');
ylabel('\theta [deg]');
grid on;
msmt_hdls = plot(1:total_num_samples,zeta*180/pi,'DisplayName','Measurements');
for ii=1:numel(msmt_hdls)
    msmt_hdls(ii).Color = [msmt_hdls(ii).Color .5];
end
utils.excludeFromLegend(msmt_hdls(2:end));
avg_hdls=plot(1:total_num_samples,zeta_avg*180/pi,'k','LineWidth',2,'DisplayName','Running Average');
utils.excludeFromLegend(avg_hdls(2:end));
curr_msmt_markers = scatter(ones(2,1),zeta(:,1),'ko','filled');
utils.excludeFromLegend(curr_msmt_markers);
curr_avg_markers = scatter(ones(2,1),zeta_avg(:,1),'ko','filled');
utils.excludeFromLegend(curr_avg_markers);
legend('Location','NorthWest');
xlim([1,total_num_samples]);
grid on;
set(gca,'xscale','log');

% Plot the solutions
subplot(212);
scatter(x_aoa(1,:)/1e3, x_aoa(2,:)/1e3, 'k^','filled','DisplayName','Sensors');
hold on;

num_tail = 25;
alpha_vec = linspace(0,1,num_tail);

hdl_est = scatter(x_est_full(1,1:num_tail)/1e3,...
                    x_est_full(2,1:num_tail)/1e3,'blue','filled',...
                    'AlphaData',alpha_vec,...
                    'MarkerFaceAlpha','flat',...
                    'DisplayName','Estimated Position');
scatter(x_tgt(1)/1e3,x_tgt(2)/1e3,'rs','filled','DisplayName','Target');
    legend('Location','NorthWest');
    xlim([-3,3]);
    ylim([0,5]);
    grid on;

this_ell = crlb_ellipses{num_tail};
hdl_ell = plot(this_ell(1,:)/1e3, this_ell(2,:)/1e3,'-.','DisplayName','CRLB');

%% Plot Frame by Frame
samples_per_frame = 1; % new frame every x samples
num_frames = 1+ceil((total_num_samples-num_tail)/samples_per_frame);

% zoom_specs = struct('K',{1,200,5000,1000,2000},...
%                      'x',{[-5,5],[-1,4], [1,3],[1.5,2.5],[1.8,2.2]},...
%                      'y',{[0,5],[1,5],[3,5],[3.6,4.4],[3.8,4.2]});
zoom_specs = struct('K',{1,30,60,100,200,400},...
                     'x',{[-5,5],[1,3],[1.5,2.5],[1.8,2.2],[1.9,2.1],[1.95,2.05]},...
                     'y',{[0,5],[3,5],[3.6,4.4],[3.8,4.2],[3.9,4.1],[3.95,4.05]});

for idx_frame = 1:num_frames
    % Find the current sample
    K = min(total_num_samples, num_tail+samples_per_frame*(idx_frame-1));

    % Update zoom
    zoom_idx = find(K > [zoom_specs(:).K],1,'last');
    subplot(212);xlim(zoom_specs(zoom_idx).x);ylim(zoom_specs(zoom_idx).y);

    % Update markers
    curr_msmt_markers.XData = K*ones(n_aoa,1);
    curr_msmt_markers.YData = zeta(:,K)*180/pi;
    curr_avg_markers.XData = K*ones(n_aoa,1);
    curr_avg_markers.YData = zeta_avg(:,K)*180/pi;

    % Update Measurements
    this_x = x_est_full(:,K+(-num_tail+1:0));

    hdl_est.XData = this_x(1,:)/1e3;
    hdl_est.YData = this_x(2,:)/1e3;

    % Update error ellipse
    this_ell = crlb_ellipses{K};

    hdl_ell.XData = this_ell(1,:)/1e3;
    hdl_ell.YData = this_ell(2,:)/1e3;

    % Redraw
    drawnow;
    pause(.01);

end