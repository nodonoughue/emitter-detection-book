function exportPlot(fig,fnm,raster_output)
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

if nargin < 3 || isempty(raster_output)
    raster_output = false;
end

if raster_output
    exportgraphics(fig,[fnm '.svg'],'ContentType','image')
else
    exportgraphics(fig,[fnm '.svg'],'ContentType', 'vector');
end