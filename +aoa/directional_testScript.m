% Example evalution of an Adcock and Rectangular-aperture DF receiver
%
% 25 July 2019
% Nicholas O'Donoughue

% Create the antenna pattern generating function
% --- NOTE --- g,g_dot take radian inputs (psi, not theta)
d_lam = .25;
[g,g_dot] = aoa.make_gain_functions('adcock',d_lam,0);
    
% Generate the angular samples and true gain values
th_true = 5;
psi_true = th_true*pi/180;
psi_res = .001; % desired resolution from multi-stage search;
                % see aoa.directional_df for details.
N = 10;
th = linspace(-180,180-360/N,N); % evenly spaced across unit circle
psi = th(:)*pi/180;
x = g((psi-psi_true)); % Actual gain values

% Set up the parameter sweep
M_vec = [1,10,100];         % Number of temporal samples at each antenna 
                            % test point
snr_db_vec = -20:2:20;      % signal to noise ratio
num_mc = 1000;              % number of monte carlo trials at each 
                            % parameter setting

% Set up output variables
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Adcock Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc / M;
    fprintf('\t M=%d',M);

    % Generate Monte Carlo Noise with unit power
    noise_base = randn(N,M,this_num_mc);
    
    % Loop over SNR levels
    for idx_snr = 1:numel(snr_db_vec)
        fprintf('.');
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        n = sqrt(noise_pwr)*noise_base;
        y = x+n;
        
        % Estimate angle of arrivel for each Monte Carlo trial
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.directional_df(y(:,:,idx_mc),psi,g,psi_res,min(psi),max(psi));
        end
        
        % Compute the RMS Error
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % Compute the CRLB for RMS Error
        crlb_psi(idx_M,idx_snr) = abs(aoa.directional_crlb(snr_db,M,g,g_dot,psi,psi_true));
    end
fprintf('done.\n');
end


fig=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Adcock DF Performance');
legend('Location','SouthWest');

text(5,17,'M=1','FontSize',10);
text(5,4,'M=10','FontSize',10);
text(5,1.4,'M=100','FontSize',10);



%% Reflector/Array Test Script
% Create the antenna pattern generating function
D_lam = 5;
[g,g_dot] = aoa.make_gain_functions('rectangular',D_lam,0);
    
% Generate the angular samples and true gain values
th_true = 5;
psi_true = th_true*pi/180;
N = 36; % number of samples
th = linspace(-180,180-360/N,N);
psi = th(:)*pi/180;
psi_res = .001; % desired resolution from multi-stage search;
                % see aoa.directional_df for details.
x = g((psi-psi_true)); % Actual gain values

% Set up the parameter sweep
M_vec = [1,10,100];     % Number of temporal samples at each antenna test 
                        % point
snr_db_vec = -10:20;    % signal to noise ratio
num_mc = 1000;          % number of monte carlo trials at each parameter 
                        % setting

% Set up output variables
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Rectangular Array Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc / M;
    fprintf('\t M=%d',M);
    
    % Generate Monte Carlo Noise with unit power
    noise_base = 1/sqrt(2)*(randn(N,M,this_num_mc)+1i*randn(N,M,this_num_mc));

    % Loop over SNR levels
    for idx_snr = 1:numel(snr_db_vec)
        fprintf('.');
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        n = sqrt(noise_pwr)*noise_base;
        y = x+n;
        
        % Estimate angle of arrival for each Monte Carlo trial
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.directional_df(y(:,:,idx_mc),psi,g,psi_res,-pi,pi);
        end
        
        % Compute the RMS Error
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % Compute the CRLB for RMS Error
        crlb_psi(idx_M,idx_snr) = abs(aoa.directional_crlb(snr_db,M,g,g_dot,psi,psi_true));
    end
fprintf('done.\n');
end


fig=figure;
for idx_M=1:numel(M_vec)
    set(gca,'ColorOrderIndex',idx_M);
    h1=semilogy(snr_db_vec,sqrt(crlb_psi(idx_M,:))*180/pi,'DisplayName','CRLB');%sprintf('CRLB, M=%d',M_vec(idx_M)));
    hold on;
    set(gca,'ColorOrderIndex',idx_M);
    h2=semilogy(snr_db_vec,rmse_psi(idx_M,:)*180/pi,'--','DisplayName','Simulation Result');%sprintf('Simulation Result, M=%d',M_vec(idx_M)));
    if idx_M~=1
        utils.excludeFromLegend(h1);
        utils.excludeFromLegend(h2);
    end
end
xlabel('$\xi$ [dB]');
ylabel('RMSE [deg]');
title('Rectangular Array DF Performance');
legend('Location','SouthWest');

text(5,17,'M=1','FontSize',10);
text(5,5,'M=10','FontSize',10);
text(5,1.7,'M=100','FontSize',10);