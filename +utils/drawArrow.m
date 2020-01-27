function hdl = drawArrow(x,y,maxHeadSize)
% hdl = drawArrow(x,y,maxHeadSize)
%
% Draw an arrow on the current axes, with the provided coordinates.  Used
% as an alternative to the annotation function.  Uses the quiver function
% to draw the arrows.
%
% INPUTS:
%   x               2-element vector of x coordinate start/end positions
%   y               2-element vector of y coordinate start/end positions
%   maxHeadSize     Maximum head size, in points
%
% OUTPUTS:
%   hdl             Handle object for arrow drawn
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(maxHeadSize)
    maxHeadSize=5;
end

hdl = quiver(x(1),     y(1),...
            x(2)-x(1),y(2)-y(1),...
            0,'MaxHeadSize',maxHeadSize );
hdl.LineStyle='-';
hdl.Color=[0 0 0];

%% Remove the arrow from legend entries
utils.excludeFromLegend(hdl);