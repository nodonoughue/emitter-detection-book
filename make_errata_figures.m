% This script contains figures that appear only in the errata
utils.initPlotSettings;

%% Small-Angle Approximation Error

d_lam = 0:.01:.5;

% Worst case error, occurs when cos(psi)= +/- 1
norm_err_sin = abs(sin(pi*d_lam) - pi*d_lam)./sin(pi*d_lam);
norm_err_cos = abs(cos(pi*d_lam) - (1-.5*(pi*d_lam).^2))./cos(pi*d_lam);

figure;
plot(d_lam,norm_err_sin,d_lam,norm_err_cos);
xlabel('$d/\lambda$');
ylabel('Normalized Error [\%]');
legend('sin(x) = x','cos(x)= 1 - x^2/2');
ylim([0 .5])
grid on;
utils.exportPlot(gcf,['figures', filesep, 'small_angle_approx']);

% Reset plot settings
utils.resetPlotSettings;