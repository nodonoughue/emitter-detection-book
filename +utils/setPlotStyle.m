function setPlotStyle(ax,styleName)
% setPlotStyle(ax,styleName)
%
% Preset script for various plot styles, to handle quick setting of
% common parameters, such as "clean" for drawings with no need to mark the
% axes.
%
% If styleName is a cell array of strings, each string is interpreted, in
% order.  This gives preference to the last element in the array, which may
% override settings called by prior elements.
%
% Where possible, the visibility of an axis element is adjusted, rather
% than its presence; to allow for toggling by later commands.
%
% Inputs:
%
%   ax          Axis handle to modify
%   styleName   String or cell array of strings.  Valid options include:
%                   'clean'     strips axes
%                   'notick'    removes numbers and ticks, but leaves axes
%                   'equal'     sets scale equal for x and y axes
%                   'widescreen'sets the plot dimensions to a 16:9 format
%                   'tight'     minimizes white space outside axes
%
% Nicholas O'Donoughue
% 1 July 2019


% Iterate over comma-separated lists; in order.
if iscell(styleName)
    cellfun(@(x) utils.setPlotStyle(ax,x),styleName);
    return;
end



switch styleName
    case 'clean'
        % Remove the tick marks and grid
%        grid(ax,'off');
%        set(ax,'xtick',[]);
%        set(ax,'ytick',[]);
%        set(ax,'xticklabel',[]);
%        set(ax,'yticklabel',[]);
        ax.Visible ='off';
        
        if ~isempty(ax.Legend)
            ax.Legend.Box = 'off';
        end
    
    case 'notick'
        set(ax,'xticklabel',[]);
        set(ax,'yticklabel',[]);
        
    case 'box only'
        set(ax,'xticklabel',[]);
        set(ax,'yticklabel',[]);
        grid(ax,'off');
        
    case 'equal'
        axis(ax,'equal');
        
        % Set figure size to square
        %posVec = get(ax.Parent,'Position');
        %posVec(4)=posVec(3);
        %set(ax.Parent,'Position',posVec);
        
        % Set axes position to square
        %posVec = get(ax,'OuterPosition');
        %posVec(4) = posVec(3);
        %set(ax,'OuterPosition',posVec);
    case 'widescreen'
        % 16:9 format
        
        % Set figure size
        posVec = get(ax.Parent,'Position');
        posVec(4) = 9*posVec(3)/16;
        set(ax.Parent,'Position',posVec);
        
        % Set axis position
        posVec = get(ax,'OuterPosition');
        posVec(4) = 9*posVec(3)/16;
        set(ax,'OuterPosition',posVec);
    case 'tight'
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
end
