function fig=ex2_2()
% fig=ex2_2()
%
% Executes Example 2.2.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Set up PFA and SNR vectors
pfa = linspace(1e-9,1,1001);
d2_vec = [1e-6,1,5,10]; % Use 1e-6 instead of 0, to avoid NaN

% Use ndgrid to simplify the code
[PFA,D2] = ndgrid(pfa,d2_vec);

% Compute the threshold eta using MATLAB's built-in error function erf(x) 
% and inverse error function erfinv(x).
eta = sqrt(2*D2).*erfinv(1-2*PFA)-D2/2;

% Compute the probability of detection
PD = .5*(1-erf((eta-D2/2)./sqrt(2*D2)));

% Plot the ROC curve
fig=figure;
plot(pfa,PD);

% Axes Labels
ylabel('$P_D$');xlabel('$P_{FA}$');

% Legend
legendStr = arrayfun(@(x) sprintf('d^2 = %d',round(x)),d2_vec,'UniformOutput',false);
legend(legendStr,'Location','SouthEast');

% Align the axes
utils.setPlotStyle(gca,{'widescreen','tight'});
xlim([0,1]);ylim([0 1]);
