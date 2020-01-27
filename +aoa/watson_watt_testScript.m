% Test script that demonstrates how to analyze a Watson Watt DF receiver
%
% 25 July 2019
% Nicholas O'Donoughue

% Generate the Signals
th_true = 45;
psi_true = th_true*pi/180;
f = 1e9;
t_samp = 1/(3*f); % ensure the Nyquist criteria is satisfied    

% Set up the parameter sweep
M_vec = [1,10,100];     % Number of temporal samples at each antenna 
                        % test point
snr_db_vec = -10:.2:20; % signal to noise ratio
num_mc = 10000;         % number of monte carlo trials at each parameter 
                        % setting

% Set up output variables
rmse_psi = zeros(numel(M_vec),numel(snr_db_vec));
crlb_psi = zeros(numel(M_vec),numel(snr_db_vec));

% Loop over parameters
fprintf('Executing Watson Watt Monte Carlo sweep...\r\n');
for idx_M = 1:numel(M_vec)
    M = M_vec(idx_M);
    this_num_mc = num_mc/M;
    fprintf('\t M=%d',M);
    
    % Generate signal vectors
    t_vec = (0:M-1)*t_samp;
    r0 = cos(2*pi*f*t_vec).';
    y0 = sind(th_true)*r0;
    x0 = cosd(th_true)*r0;
    
    % Generate Monte Carlo Noise with unit power
    noise_base_r = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);
    noise_base_x = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);
    noise_base_y = sqrt(sum(abs(r0).^2)/M)*randn(M,this_num_mc);

    % Loop over SNR levels
    for idx_snr = 1:numel(snr_db_vec)
        if mod(idx_snr,10)==0
            fprintf('.');
        end
        
        % Compute noise power, scale base noise
        snr_db = snr_db_vec(idx_snr);
        noise_pwr = 10.^(-snr_db/10);
        
        % Generate noisy measurement
        r = r0 + noise_base_r*sqrt(noise_pwr);
        y = y0 + noise_base_y*sqrt(noise_pwr);
        x = x0 + noise_base_x*sqrt(noise_pwr);
        
        % Compute the estimate for each Monte Carlo trial
        psi_est = zeros(1,this_num_mc);
        for idx_mc = 1:this_num_mc
            psi_est(idx_mc) = aoa.watson_watt_df(r(:,idx_mc),x(:,idx_mc),y(:,idx_mc));
        end

        % Compute RMS Error
        rmse_psi(idx_M,idx_snr) = sqrt(sum(abs((psi_est-psi_true)).^2)/this_num_mc);
        
        % Compute CRLB for RMS Error
        crlb_psi(idx_M,idx_snr) = abs(aoa.watson_watt_crlb(snr_db,M));
    end
fprintf('done.\n');
end


figure;
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
title('Watson Watt DF Performance');
legend('Location','SouthWest');

text(5,17,'M=1','FontSize',10);
text(5,5,'M=10','FontSize',10);
text(5,1.7,'M=100','FontSize',10);