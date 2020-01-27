function excludeFromLegend(graphObj)
% excludeFromLegend(graphObj)
%
% Modifies the supplied graphics object such that it will not appead in a
% legend; useful for plotting annotations that do not belong in the legend,
% such as dashed lines between a label and the appropriate line or
% marker.
%
% The handles are returned from the plotting command that was used to
% generate them, such as
%   h = plot(...)
% or
%   h = text(...)
%
% Alternatively, an axis object can be searched with the findall or findobj
% commands, such as:
%   h1 = findall(gcf);
%   h2 = findall(gcf,'Type','text');
%   h3 = findobj(gcf,'Type','text');
%
% If multiple graph objects are input the command is applied to each.
%
% Inputs:
%   graphObj        Graphics object handles (scalar or array)
%
% Nicholas O'Donoughue
% 1 July 2019

if numel(graphObj)>1
    arrayfun(@(x) utils.excludeFromLegend(x),graphObj);
    return;
end

% Graphics Objects have an Annotation property.
% Within the Annotation object, there is a LegendInformation Property.
% The LegendInformation object has an IconDisplayStyle property.
%
% Set that IconDisplayProperty to 'off' and a graphics object will be
% excluded from legends
graphObj.Annotation.LegendInformation.IconDisplayStyle = 'off';