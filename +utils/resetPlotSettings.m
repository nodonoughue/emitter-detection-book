% resetPlotSettings
%
% Resets MATLAB's graphics defaults, to undo all of the changes made by
% initPlotSettings.
%
% Nicholas O'Donoughue
% 11 March 2020

% Serif font family
set(groot,'DefaultAxesFontName','default');
set(groot,'DefaultAxesFontSize','default');
set(groot,'DefaultTextFontSize','default');

% LaTeX Interpretation
set(groot,'DefaultTextInterpreter','default');

% Title Size & Weight
set(groot,'DefaultAxesTitleFontSizeMultiplier','default');
set(groot,'DefaultAxesTitleFontWeight','default');
set(groot,'DefaultAxesLabelFontSizeMultiplier','default');

% Line Width, Color, and Style
set(groot,'DefaultAxesColorOrder','default');
set(groot,'DefaultAxesLineStyleOrder','default');
set(groot,'DefaultLineLineWidth','default');

% Background Color
set(groot,'DefaultAxesColor','default');

% Grid
set(groot,'DefaultAxesXGrid','default','DefaultAxesYGrid','default');
set(groot,'DefaultAxesGridColor','default');
set(groot,'DefaultAxesMinorGridColor','default');
set(groot,'DefaultAxesGridAlpha','default');
set(groot,'DefaultAxesLineWidth','default');

% Figure Size
set(groot,'DefaultAxesUnits','default');
set(groot,'DefaultAxesOuterPosition','default');
set(groot,'DefaultFigurePosition','default');
set(groot,'DefaultFigureUnits','default');
%set(groot,'DefaultFigureOuterPosition','default');
