function fig=ex5_3()
% fig=ex5_3()
%
% Executes Example 5.3.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Define input parameters
Bs = [1e5:1e5:1e7];
Bhop = 4e9;
Thop = [1e-3,1e-4,1e-5]';
duty = .2;
tp = duty*Thop;

Tdwell = 1./Bs;
Nscan = tp./Tdwell;
Br = max(Bs,Bhop./Nscan);

fig=figure;
colors = get(0,'DefaultAxesColorOrder');
set(0,'DefaultAxesColorOrder',colors(1,:));
loglog(Bs/1e6,Br/1e6);
xlabel('Frequency Resolution ($\delta_f$) [MHz]');
ylabel('Receiver Bandwidth ($B_r$) [MHz]');
legend(arrayfun(@(x) sprintf('t_p=%.0f \\mu s',x*1e6),tp,'UniformOutput',false));

utils.setPlotStyle(gca,{'widescreen','tight'});