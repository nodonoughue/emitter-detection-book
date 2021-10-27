% Define source positions
r_min = 1e3;
r_max = 200e3;
M = 2001;
yvec = linspace(r_min,r_max,M);
x_source = [zeros(size(yvec));yvec];

% Define Sensor Positions
for baseline = [10e3,30e3]
    
    nSensors = 2;
    x_sensor = [-baseline/2, baseline/2;0 0];
    
    % Define Sensor Performance
    c_psi_vec = .1:.05:5;
    
    % Initialize output
    cep50 = zeros(numel(c_psi_vec),M);
    
    for idx = 1:numel(c_psi_vec)
        c_psi = c_psi_vec(idx)*pi/180;
        Caoa = c_psi^2 * eye(nSensors);

        % Compute CRLB
        warning('off','MATLAB:nearlySingularMatrix'); % We know the problem is ill-defined, deactivate the warning
        crlb = triang.computeCRLB(x_sensor,x_source,Caoa); % Ndim x Ndim x M^2
        cep50(idx,:) = utils.computeCEP50(crlb);
        warning('on','MATLAB:nearlySingularMatrix'); % Reactivate the singular matrix warning
    end
    
    figure;
    levels=[0:1:10,20:10:100];
    [C,h]=contourf(yvec/1e3,c_psi_vec,cep50/1e3,levels);
    clabel(C,h);
    colormap(flipud(viridis(100)));
    caxis([0 50])
    h=colorbar;
    ylabel(h,'CEP_{50} [km]')
    xlabel('Distance from Base [km]');
    ylabel('Sensor DF Accuracy [deg]');
    grid on;
    title(sprintf('Baseline = %d km',baseline/1e3))
end