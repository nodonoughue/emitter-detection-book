function fig = book2_ex6_1()
% fig=book2_ex6_1()
%
% Executes Example 6.1 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 28 November 2021

%% Set up sensors
x_aoa = [2, 2, 0; 2, -1, 0];
n_aoa = size(x_aoa,2);

x_tgt = [5; 3];

% Define received signals and covariance matrix
alpha = [5, 10, -5]*pi/180; % AOA bias
psi = triang.measurement(x_aoa, x_tgt);
psi_bias = triang.measurement(x_aoa, x_tgt, false, alpha);



%% Plot the scenario
fig = figure;
scatter(x_aoa(1,:),x_aoa(2,:),'s','filled','DisplayName','Sensors')
hold on;
scatter(x_tgt(1), x_tgt(2),'^','filled','DisplayName','Target')
grid on;

% Draw the LOBs
clr_idx = 1; % grab the current color index
xy_lob = triang.drawLob(x_aoa, psi, x_tgt,1.5);
for idx=1:n_aoa
    set(gca,'ColorOrderIndex',clr_idx); % reset to the same point; so the lobs have repeated colors
    hdl=plot(xy_lob(1,:,idx),xy_lob(2,:,idx),'DisplayName','LOB (w/o bias)');
    if idx > 1
        utils.excludeFromLegend(hdl);
    end
end

clr_idx=2;
xy_lob_bias = triang.drawLob(x_aoa, psi_bias, x_tgt, 1.5);
for idx=1:n_aoa
    set(gca,'ColorOrderIndex',clr_idx); % reset to the same point; so the lobs have repeated colors
    hdl=plot(xy_lob_bias(1,:,idx), xy_lob_bias(2,:,idx),'--','DisplayName','LOB (w/bias)');
    if idx > 1
        utils.excludeFromLegend(hdl);
    end
end

legend('Location','NorthWest');
xlim([-1 6]);
ylim([-1 5]);

utils.setPlotStyle(gca,'widescreen');
