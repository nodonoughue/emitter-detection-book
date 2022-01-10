function exportPlot(fig,fnm,formats)
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

%% Parse Inputs
if nargin < 3 || isempty(formats)
    formats = {'eps','png','fig'};
end

if exist('OCTAVE_VERSION','builtin')==0
  % Only do this command in MATLAB
  fig.InvertHardcopy = 'off';
end

if ~iscell(formats)
    formats = {formats};
end

for idx=1:numel(formats)
    if strcmpi(formats{idx},'fig')==1
        saveas(fig,[fnm '.fig']);
    elseif strcmpi(formats{idx},'eps')==1 || strcmpi(formats{idx},'epsc')==1
        print(fig,fnm,'-depsc','-r1200','-painters');
    elseif strcmpi(formats{idx},'png')==1
    	print(fig,fnm,'-dpng','-r1200');
    end
end
