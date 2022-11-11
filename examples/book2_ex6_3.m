function fig = book2_ex6_3()
% fig=book2_ex6_3()
%
% Executes Example 6.3 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 8 February 2022

%% Set up sensors
x_aoa = [2, 2, 0; 2,  -1, 0];
[num_dims,n_aoa] = size(x_aoa);

x_tdoa = [0, 2, 0; 2, -1, 0]; % avg position (reported)
n_tdoa = size(x_tdoa,2);

% Find matching sensors
dist = utils.rng(x_aoa,x_tdoa);
idx_pairs = find(dist(:)==0);
[idx_aoa,idx_tdoa] = ind2sub(size(dist),idx_pairs);

% Build position covariance matrix
C_pos_1d = eye(n_tdoa+n_aoa); % position covar; all are IID
for ii = 1:numel(idx_pairs)
    this_idx_aoa = idx_aoa(ii);
    this_idx_tdoa = idx_tdoa(ii);
    C_pos_1d(this_idx_aoa,n_aoa+this_idx_tdoa) = 1;
    C_pos_1d(n_aoa+this_idx_tdoa,this_idx_aoa) = 1;
end
% At this point, C_pos_1d is the correlation coefficient of all
% sensor position errors in a given dimension.  We must now expand
% the matrix to account for (a) actual position variance and (b) the
% number of spatial dimensions in the problem.
C_beta_single = .1*eye(num_dims);
C_beta = kron(C_pos_1d,C_beta_single);

% Generate a random set of AOA and TDOA positions
% L=chol(C_beta,'lower'); -- This will fail, because C_beta is not positive
% definite (it has some eigenvalues that are zero)
[U,S,~] = svd(C_beta);
epsilon = reshape(U*sqrt(S)*randn(num_dims*(n_aoa+n_tdoa),1),num_dims,n_aoa+n_tdoa);

% Grab the position offsets and add to the reported TDOA and AOA positions
beta_aoa = x_aoa + epsilon(:,1:n_aoa);
beta_tdoa = x_tdoa + epsilon(:,n_aoa + (1:n_tdoa));

% Let's verify that sensors 2 and 4 are still colocated
dist_perturbed = utils.rng(beta_aoa,beta_tdoa);
assert(all(dist_perturbed(idx_pairs)<1e-6),'Error generating correlated sensor perturbations.');

%% Generate Measurements
x_tgt = [6; 3];

alpha_aoa = [5, 10, -5]*pi/180; % AOA bias
zeta = hybrid.measurement(x_aoa,x_tdoa,[],[], x_tgt, n_tdoa); % free of pos unc and bias
zeta_unc = hybrid.measurement(beta_aoa,beta_tdoa,[],[],x_tgt,n_tdoa); % with pos unc, no bias
zeta_unc_bias = hybrid.measurement(beta_aoa,beta_tdoa,[],[], x_tgt, n_tdoa,[],[],alpha_aoa); % with pos unc and bias

fprintf('Measurements from ideal sensors (AOA, AOA, AOA, RDOA, RDOA):\n[')
fprintf('%.2f, ',zeta);
fprintf(']\nWith pos unc:\n[');
fprintf('%.2f, ',zeta_unc);
fprintf(']\nWith pos unc and bias:\n[');
fprintf('%.2f, ',zeta_unc_bias);
fprintf(']\n');

%% Plot Scenario
fig = figure;
scatter([x_aoa(1,:),x_tdoa(1,:)],[x_aoa(2,:),x_tdoa(2,:)],'s','filled','DisplayName','Sensors (nominal positions)')
hold on;
scatter([beta_aoa(1,:), beta_tdoa(1,:)],[beta_aoa(2,:),beta_tdoa(2,:)],'o','filled','DisplayName','Sensors (true positions)')
scatter(x_tgt(1), x_tgt(2),'^','filled','DisplayName','Target')
grid on;

% Draw the Isochrones and LOBs -- Truth
xy_lob = triang.drawLob(x_aoa, zeta_unc_bias(1:n_aoa), x_tgt, 1.5);
% xy_lob_bias = triang.drawLob(beta_aoa,zeta_unc_bias(1:n_aoa), x_tgt, 1.5);
for idx=1:n_aoa
    set(gca,'ColorOrderIndex',3); % reset to the same point; so the lobs have repeated colors
    hdl=plot(xy_lob(1,:,idx),xy_lob(2,:,idx),'DisplayName','LOB (nominal positions)');

%     set(gca,'ColorOrderIndex',idx); % reset to the same point; so the lobs have repeated colors
%     hdl_bias=plot(xy_lob_bias(1,:,idx),xy_lob_bias(2,:,idx),'-.','DisplayName','Isochrone/LOB (w/pos.unc. and bias)');
    if idx > 1
        utils.excludeFromLegend(hdl);
    end
end
for idx=1:n_tdoa-1
    xy_iso = tdoa.drawIsochrone(x_tdoa(:,end),x_tdoa(:,idx),zeta_unc_bias(n_aoa+idx),101,8);
%     xy_iso_bias= tdoa.drawIsochrone(beta_tdoa(:,end), beta_tdoa(:,idx),zeta_unc_bias(n_aoa+idx),101,8);

    set(gca,'ColorOrderIndex',4); % reset to the same point; so the lobs have repeated colors
    hdl=plot(xy_iso(1,:),xy_iso(2,:),'DisplayName','Isochrone (nominal positions)');
%     set(gca,'ColorOrderIndex',n_aoa+idx); % reset to the same point; so the lobs have repeated colors
%     hdl_bias=plot(xy_iso_bias(1,:),xy_iso_bias(2,:),'-.','DisplayName','Isochrone (w/pos. unc. and bias)');
    if idx > 1
        utils.excludeFromLegend(hdl);
    end
end

legend('Location','SouthEast');
xlim([-1 8]);
ylim([-3 4]);

utils.setPlotStyle(gca,'widescreen');