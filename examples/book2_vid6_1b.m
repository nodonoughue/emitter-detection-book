doFullColor = true;
utils.initPlotSettings;
colors=get(groot,'DefaultAxesColorOrder');

%% Figures 6.3a and 6.3b, Impact of Sensor Position Errors

x_tdoa = [-1, 1, 3;0, 0, .5];
x_target = [0; 5];

%% Compute CRLB
[n_dim, n_tdoa] = size(x_tdoa);
pos_var = 1;
C_pos = pos_var*eye(n_dim*n_tdoa);

x_vec = linspace(-6,6,101);
y_vec = linspace(0,11,101);
[XX,YY] = meshgrid(x_vec,y_vec);
x_grid = [XX(:),YY(:)]';

n_grid = size(x_grid,2);
crlb = zeros(n_dim,n_dim,n_grid);

J = tdoa.grad_b(x_tdoa,x_grid);

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
set(gca,'ydir','normal');
hold on;
scatter(x_tdoa(1,:),x_tdoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');


%% X Error
C_single_sensor = diag([pos_var,pos_var/1000]);
C_x = kron(eye(n_tdoa),C_single_sensor);

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
set(gca,'ydir','normal');
hold on;
scatter(x_tdoa(1,:),x_tdoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');

%% Y Error
C_single_sensor = diag([pos_var/1000,pos_var]);
C_y = kron(eye(n_tdoa),C_single_sensor);

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
set(gca,'ydir','normal');
hold on;
scatter(x_tdoa(1,:),x_tdoa(2,:),'v','filled','DisplayName','Sensors');
scatter(x_target(1),x_target(2),'^','filled','DisplayName','Target');
grid on;
legend('Location','NorthEast');