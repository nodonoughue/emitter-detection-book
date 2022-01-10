% Test ellipse approximation
close all;

rad_max = 5e3; % max radius
ecc_vec = .1:.1:10; % ratio of min to max eigenvalue
cep_vec = 10.^((-10:.5:10)/10)*rad_max; % vector of cep50 for each test

n_mc = 100000;
n_gamma_pts = 101;

n_ecc = numel(ecc_vec);
n_cep = numel(cep_vec);

pd_mc = zeros(n_ecc, n_cep);
pd_approx = zeros(n_ecc, n_cep);close all


for i_ecc = 1:n_ecc
    ecc = ecc_vec(i_ecc);
    
    % Construct a 2D covariance matrix
    cov = diag([1, ecc]);
    
    for i_cep = 1:n_cep
        cep = cep_vec(i_cep);
        
        % Scale the covariance matrix
        scale = cep / utils.computeCEP50(cov);
        cov = scale^2 * cov;
        
        % Monte Carlo Trial
        lower = chol(cov,'lower');
        x = lower*randn(2,n_mc);
        
        % Compute radius and compare to desired
        radius = sqrt(sum(abs(x).^2,1));
        count = sum(radius < rad_max);
        pd_mc(i_ecc, i_cep) = count/n_mc;
        
        % Use approximation method
%         lam = eig(cov);
%         a2 = lam(1); b2 = lam(2);
%         
% %         th_vec = linspace(0,pi/2,n_gamma_pts);
% %         rad_sq_vec = a2*b2 ./ (a2*sin(th_vec).^2 + b2*cos(th_vec).^2); 
%         gamma_sqrt = rad_max ./ sqrt(lam);
%         gamma_vec = linspace(gamma_sqrt(1),gamma_sqrt(2),n_gamma_pts).^2;
%         prob_handoff_vec = utils.computeRMSEConfInterval(gamma_vec);
%         pd_approx(i_ecc, i_cep) = mean(prob_handoff_vec);
        pd_approx(i_ecc, i_cep) = utils.ellipseCDF(cov, rad_max);
    end
    
%     figure;
%     scatter(pd_mc(i_ecc,:),pd_approx(i_ecc,:));
%     hold on;plot([0,1],[0,1],'k-.');
%     xlabel('Monte Carlo');
%     ylabel('Approximation');
%     title(sprintf('Eccentricity: %.2f',ecc));
%     grid on;
end

figure;
scatter(pd_mc(:),pd_approx(:));
hold on;plot([0, 1],[0, 1],'k-.');
xlabel('Monte Carlo');
ylabel('Approximation');
grid on;


rmse_cep = sqrt(sum((pd_mc - pd_approx).^2,1));
rmse_ecc = sqrt(sum((pd_mc - pd_approx).^2,2));
rmse_full = sqrt(sum((pd_mc(:) - pd_approx(:)).^2));

figure;
plot(cep_vec/rad_max,rmse_cep);grid on;
ylim([0, .05]);
title('RMSE of P_D Approximation, as a function of CEP Ratio');
xlabel('CEP/R');
ylabel('RMSE');

figure;
plot(ecc_vec,rmse_ecc);grid on;
ylim([0, .05]);
title('RMSE of P_D Approximation, as a function of Eccentricity');
xlabel('Covariance Eccentricity (\sigma_1^2/\sigma_2^2)');
ylabel('RMSE');

