function exportPlot(fig,fnm)
% exportPlot(fig,fnm)
%
% Exports a plot to file using the desired resolution and format for
% inclusion in a textbook.
%
% Currently exports both .Fig and .EPS files with 1200 dpi resolution.
%
% Inputs:
%   fig         Figure handle to export
%   fnm         Filename, including directory
%
% Nicholas O'Donoughue
% 1 July 2019

fig.InvertHardcopy = 'off';

saveas(fig,[fnm '.fig']);
print(fig,fnm,'-depsc','-r1200','-painters');
print(fig,fnm,'-dpng','-r1200');
