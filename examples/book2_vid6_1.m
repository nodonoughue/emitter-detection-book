doFullColor = true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

%% Figures 6.3a and 6.3b, Impact of Sensor Position Errors

x_aoa = [-1, 1;0, 0];
x_target = [0; 5];
x_err = [1;0];
x_err_offaxis = [0;1];

x_aoa_un = x_aoa + x_err;
x_aoa_nonun = x_aoa + x_err*[-1 1];
x_aoa_offaxis = x_aoa + x_err_offaxis*[-1 1];

% Solve True LOBs and AOAs
lob = x_target - x_aoa;
psi = triang.measurement(x_aoa, x_target);
r = utils.rng(x_aoa, x_target);

% Solve Erroneous AOAs
x_tgt_un = utils.find_intersect(x_aoa_un(:,1),psi(1),x_aoa_un(:,2),psi(2));
x_tgt_nonun = utils.find_intersect(x_aoa_nonun(:,1),psi(1),x_aoa_nonun(:,2),psi(2));
x_tgt_offaxis = utils.find_intersect(x_aoa_offaxis(:,1),psi(1),x_aoa_offaxis(:,2),psi(2));

% num_dim x num_tx x 2
lob_zero = [cos(psi), sin(psi)]'.*reshape([zeros(2,1),3*r],1,2,2);

lob_true = x_aoa + lob_zero;
lob_un = x_aoa_un + lob_zero;
lob_nonun = x_aoa_nonun + lob_zero;
lob_offaxis = x_aoa_offaxis + lob_zero;

% Draw Uniform Offset case
fig3a=figure;
hdls = plot(squeeze(lob_true(1,:,:))',...
            squeeze(lob_true(2,:,:))','-','Color',colors(1,:),...
            'DisplayName','True LOB');
utils.excludeFromLegend(hdls(2:end));
hold on;
hdls = plot(squeeze(lob_un(1,:,:))',...
            squeeze(lob_un(2,:,:))','--','Color',colors(2,:),...
            'DisplayName','Perceived LOB');
utils.excludeFromLegend(hdls(2:end));
plot(x_aoa(1,:),x_aoa(2,:),'^','Color',colors(1,:),'DisplayName','True Sensor Position');
plot(x_aoa_un(1,:),x_aoa_un(2,:),'v','Color',colors(2,:),'DisplayName','Est. Sensor Position');
plot(x_target(1),x_target(2),'o','Color',colors(1,:),'DisplayName','Target');
plot(x_tgt_un(1),x_tgt_un(2),'o','Color',colors(2,:),'DisplayName','Est. Target');
xlim([-6 6]);
ylim([0 11]);
legend('Location','NorthWest');

utils.setPlotStyle(gca,{'tight','clean','equal'});

% Draw Non-Uniform Offset case
fig3b=figure;
plot(squeeze(lob_true(1,:,:))',...
     squeeze(lob_true(2,:,:))','-','Color',colors(1,:),'DisplayName','True LOB');
hold on;
plot(squeeze(lob_nonun(1,:,:))',...
     squeeze(lob_nonun(2,:,:))','--','Color',colors(2,:),'DisplayName','Perceived LOB');
plot(x_aoa(1,:),x_aoa(2,:),'^','Color',colors(1,:),'DisplayName','True Sensor Position');
plot(x_aoa_nonun(1,:),x_aoa_nonun(2,:),'v','Color',colors(2,:),'DisplayName','Est. Sensor Position');
plot(x_target(1),x_target(2),'o','Color',colors(1,:),'DisplayName','Target');
plot(x_tgt_nonun(1),x_tgt_nonun(2),'o','Color',colors(2,:),'DisplayName','Est. Target');
% legend('Location','NorthWest');
xlim([-6 6]);
ylim([0 11]);

utils.setPlotStyle(gca,{'tight','clean','equal'});


% Draw Vertical Offset case
fig3c=figure;
plot(squeeze(lob_true(1,:,:))',...
     squeeze(lob_true(2,:,:))','-','Color',colors(1,:),'DisplayName','True LOB');
hold on;
plot(squeeze(lob_offaxis(1,:,:))',...
     squeeze(lob_offaxis(2,:,:))','--','Color',colors(2,:),'DisplayName','Perceived LOB');
plot(x_aoa(1,:),x_aoa(2,:),'^','Color',colors(1,:),'DisplayName','True Sensor Position');
plot(x_aoa_offaxis(1,:),x_aoa_offaxis(2,:),'v','Color',colors(2,:),'DisplayName','Est. Sensor Position');
plot(x_target(1),x_target(2),'o','Color',colors(1,:),'DisplayName','Target');
plot(x_tgt_offaxis(1),x_tgt_offaxis(2),'o','Color',colors(3,:),'DisplayName','Est. Target');
% legend('Location','NorthWest');
xlim([-6 6]);
ylim([0 11]);

utils.setPlotStyle(gca,{'tight','clean','equal'});

%% Compute CRLB
[n_dim, n_aoa] = size(x_aoa);
pos_var = 1;
C_pos = pos_var*eye(n_dim*n_aoa);

x_vec = linspace(-6,6,101);
y_vec = linspace(0,11,101);
[XX,YY] = meshgrid(x_vec,y_vec);
x_grid = [XX(:),YY(:)]';

n_grid = size(x_grid,2);
crlb = zeros(n_dim,n_dim,n_grid);

J = triang.grad_b(x_aoa,x_grid);

for idx_pos = 1:n_grid
    this_J = squeeze(J(:,:,idx_pos));
    crlb(:,:,idx_pos) = inv(this_J'/C_pos * this_J);
end

crlb_rmse = reshape(sqrt(crlb(1,1,:)+crlb(2,2,:)),size(XX));
levels=[.01,.02,.05,.1,.2,.5,1,2,5,10,20,50,100];

figure;
title(sprintf('Uniform X/Y Error; Variance=%.2f',pos_var));
%imagesc(x_vec,y_vec,crlb_rmse);
[C,h] = contourf(XX,YY,crlb_rmse,levels);
clabel(C,h);
utils.excludeFromLegend(h);
colorbar;
caxis([0,10]);
set(gca,'ydir','normal');
hold on;
scatter(x_aoa(1,:),x_aoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');


%% X Error
C_single_sensor = diag([pos_var,pos_var/1000]);
C_x = kron(eye(n_aoa),C_single_sensor);

crlb_x = zeros(n_dim,n_dim,n_grid);

for idx_pos = 1:n_grid
    this_J = squeeze(J(:,:,idx_pos));
    crlb_x(:,:,idx_pos) = inv(this_J'/C_x * this_J);
end

crlb_x_rmse = reshape(sqrt(crlb_x(1,1,:)+crlb_x(2,2,:)),size(XX));
figure;
title(sprintf('X Error; Variance=%.2f',pos_var));
%imagesc(x_vec,y_vec,crlb_rmse);
[C,h] = contourf(XX,YY,crlb_x_rmse,levels);
clabel(C,h);
utils.excludeFromLegend(h);
colorbar;
caxis([0,10]);
set(gca,'ydir','normal');
hold on;
scatter(x_aoa(1,:),x_aoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');

%% Y Error
C_single_sensor = diag([pos_var/1000,pos_var]);
C_y = kron(eye(n_aoa),C_single_sensor);

crlb_y = zeros(n_dim,n_dim,n_grid);

for idx_pos = 1:n_grid
    this_J = squeeze(J(:,:,idx_pos));
    crlb_y(:,:,idx_pos) = inv(this_J'/C_y * this_J);
end

crlb_y_rmse = reshape(sqrt(crlb_y(1,1,:)+crlb_y(2,2,:)),size(XX));
figure;
title(sprintf('Y Error; Variance=%.2f',pos_var));
%imagesc(x_vec,y_vec,crlb_rmse);
[C,h] = contourf(XX,YY,crlb_y_rmse,levels);
clabel(C,h);
utils.excludeFromLegend(h);
colorbar;
caxis([0,10]);
set(gca,'ydir','normal');
hold on;
scatter(x_aoa(1,:),x_aoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');

