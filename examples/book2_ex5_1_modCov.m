% Executes a modified Example 5.1 from Practical Geolocation for Electronic 
% Warfare with MATLAB, used in Video 5.1.
%
% Nicholas O'Donoughue
% 7 February 2022

%% Set up sensors
x_aoa = [-2, 2; 0, 0];
[~,n_aoa] = size(x_aoa);

% Define received signals and covariance matrix
psi = [80; 87]*pi/180;
C = diag([.05, .2]); % modify variance of first to be more reliable
x_init = [0;1]; % initial guess

%% Plot the scenario
fig1 = figure;
plot(x_aoa(1,:),x_aoa(2,:),'+','DisplayName','Sensors')
hold on;
set(gca,'ColorOrderIndex',1);
grid on;

% Draw the LOBs
lob_len = 35;
lob = lob_len * [cos(psi), sin(psi)]';
x_lob = x_aoa(1,:) + [zeros(1,n_aoa); lob(1,:)];
y_lob = x_aoa(2,:) + [zeros(1,n_aoa); lob(2,:)];
hdl=plot(x_lob,y_lob,'k-','DisplayName','LOBs');
utils.excludeFromLegend(hdl(2:end));
utils.setPlotStyle(gca,'widescreen');

%% Gradient Descent Solution (unconstrained)
[x_gd, x_gd_full] = triang.gdSoln(x_aoa, psi, C, x_init);

%% Gradient Descent Solution (constrained)
y_soln = 25;
[a, ~] = utils.constraints.fixedCartesian('y',y_soln);
[x_gd_const,x_gd_full_const] = triang.gdSolnFixed(x_aoa, psi, C, x_init, a);

[x_ls, x_ls_full] = triang.lsSoln(x_aoa, psi, C, x_init);
[x_ls_const,x_ls_full_const] = triang.lsSolnFixed(x_aoa, psi, C, x_init, a);


%% Print Results
fprintf('Unconstrained Solution: (%.2f, %.2f)\n', x_gd(1), x_gd(2));
fprintf('Constrained Solution:   (%.2f, %.2f)\n', x_gd_const(1), x_gd_const(2));

fprintf('Unconstrained LS Solution: (%.2f, %.2f)\n', x_ls(1), x_ls(2));
fprintf('Constrained LS Solution:   (%.2f, %.2f)\n', x_ls_const(1), x_ls_const(2));

%% Plot with Solutions
fig2 = figure;
plot(x_aoa(1,:),x_aoa(2,:),'k+','DisplayName','Sensors')
hold on;
set(gca,'ColorOrderIndex',1);
grid on;

% Draw the LOBs
lob_len = 35;
lob = lob_len * [cos(psi), sin(psi)]';
x_lob = x_aoa(1,:) + [zeros(1,n_aoa); lob(1,:)];
y_lob = x_aoa(2,:) + [zeros(1,n_aoa); lob(2,:)];
hdl=plot(x_lob,y_lob,'k-','DisplayName','LOBs');
utils.excludeFromLegend(hdl(2:end));

% Unconstrained Solution
hdl=plot(x_gd_full(1,:),x_gd_full(2,:),'--',...
     'DisplayName','GD (unconstrained)');
plot(x_gd(1),x_gd(2),'-.o','Color',hdl.Color,...
    'DisplayName','GD (unconstrained)');
utils.excludeFromLegend(hdl);

% Constrained Solution
hdl=plot(x_gd_full_const(1,:), x_gd_full_const(2,:),'--');
plot(x_gd_const(1),x_gd_const(2),'-.s','Color',hdl.Color,...
    'DisplayName','GD (constrained)');
utils.excludeFromLegend(hdl);
legend('Location','NorthWest');
utils.setPlotStyle(gca,'widescreen');

% Unconstrained Solution (LS)
hdl=plot(x_ls_full(1,:),x_ls_full(2,:),':',...
     'DisplayName','LS (unconstrained)');
plot(x_ls(1),x_ls(2),':o','Color',hdl.Color,...
    'DisplayName','LS (unconstrained)');
utils.excludeFromLegend(hdl);

% Constrained Solution (LS)
hdl=plot(x_ls_full_const(1,:), x_ls_full_const(2,:),':');
plot(x_ls_const(1),x_ls_const(2),':s','Color',hdl.Color,...
    'DisplayName','LS (constrained)');
utils.excludeFromLegend(hdl);
legend('Location','NorthWest');
utils.setPlotStyle(gca,'widescreen');

