function figs = book2_ex7_1()
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

% Define sensor positions
x_aoa = [0, 1e3;
         0, 0];
x_tgt = [2e3;
         4e3];

% Define sensor accuracy
n_sensors = size(x_aoa,2);
sigma_theta = 5;
sigma_psi = sigma_theta * pi/180;
C = sigma_psi^2 * eye(n_sensors);

% Compute CRLB
crlb = triang.computeCRLB(x_aoa, x_tgt, C);
cep50=utils.computeCEP50(crlb);
fprintf('CEP50: %.2f km\n',cep50/1e3);

cep50_desired = 100;
est_k_min = ceil((cep50/cep50_desired).^2);
fprintf('Estimate %d samples required.\n',est_k_min);

% Iterate over number of samples
num_samples = [1:1000,2000:1000:10000];
sigma_theta_vec = [5, 10, 30];

cep_vec = zeros(numel(num_samples), numel(sigma_theta_vec));
for idx_n = 1:numel(num_samples)
    this_num_samples = num_samples(idx_n);
    for idx_s = 1:numel(sigma_theta_vec)
        this_sigma_psi = sigma_theta_vec(idx_s)*pi/180;
        this_C = this_sigma_psi^2/this_num_samples*eye(n_sensors);

        this_crlb = triang.computeCRLB(x_aoa, x_tgt, this_C);
        cep_vec(idx_n, idx_s) = utils.computeCEP50(this_crlb);
    end
end

fig=figure;
for idx_s = 1:numel(sigma_theta_vec)
    loglog(num_samples, cep_vec(:,idx_s),'DisplayName',...
             sprintf('\\sigma_\\theta=%d^\\circ',sigma_theta_vec(idx_s)));
    hold on;
end
plot(num_samples, 100*ones(size(num_samples)), 'k--', 'DisplayName','CEP=100 m');
plot(est_k_min*[1,1],[1e1,1e4],'k--','DisplayName',sprintf('K=%d',est_k_min));
ylim([10, 10e3]);
xlabel('Number of Samples [K]');
ylabel('$CEP_{50}$ [m]');
grid on;
legend('Location','SouthWest');
utils.setPlotStyle(gca,{'widescreen'});

% Determine when CEP50 crosses below 100 m
desired_cep = 100;
first_good_sample = find(cep_vec(:,1) <= desired_cep, 1, 'first');
if isempty(first_good_sample)
    fprintf('More than %d samples required to achieve %.2f m CEP50.\n', max(num_samples), desired_cep);
else
    fprintf('%d samples required to achieve %.2f m CEP50.\n', num_samples(first_good_sample), desired_cep);
end

% Plot for AOR
n_pts = 101;
x_vec = linspace(-5e3,5e3,n_pts);
y_vec = linspace(1e3,5e3,n_pts);
[xx,yy] = meshgrid(x_vec,y_vec);
x_aor = [xx(:), yy(:)]';

aor_dims = size(xx);
% n_aor = size(x_aor,2);

k_vec = [10, 100];
cmap_lim = [0, 5];

clear figs;
figs(numel(k_vec)+1) = fig;
figs(1) = fig;
for idx_k = 1:numel(k_vec)
    this_k = k_vec(idx_k);

    this_crlb = triang.computeCRLB(x_aoa, x_aor, C/this_k);
    cep = reshape(utils.computeCEP50(this_crlb), aor_dims);

    this_fig = figure;
    figs(idx_k+1) = this_fig;
    imagesc(x_vec/1e3, y_vec/1e3, cep/1e3);
    hold on;
    [Ch,h] = contour(x_vec/1e3, y_vec/1e3, cep/1e3,'k-.');
    clabel(Ch,h);
    xlabel('x [km]');
    ylabel('y [km]');
    colorbar;
    colormap(flipud(utils.viridis));
    caxis(cmap_lim);
    plot(x_aoa(1,:)/1e3, x_aoa(2,:)/1e3,'k^','DisplayName','Sensors');
    set(gca,'ydir','normal');
    ylim([-.5 5.5]);

    utils.setPlotStyle(gca,{'equal'});
end


