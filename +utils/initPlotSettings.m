% initPlotSettings
%
% Initializes plotting preferences, using global values, to a standard
% that is appropriate for publication in a textbook.
%
% Based loosely on the Seaborn-paper plotting style in the Python package
% matplotlib.
%
% Nicholas O'Donoughue
% 1 July 2019

if ~exist('doFullColor','var')
    % By default, generate grayscale images.  But, if doFullColor is
    % defined, and set to true, we'll generate full color images.
    doFullColor=false;
end

% Serif font family
set(groot,'DefaultAxesFontName','Times New Roman');
set(groot,'DefaultAxesFontSize',10);
set(groot,'DefaultTextFontSize',12);

% LaTeX Interpretation
set(groot,'DefaultTextInterpreter','latex'); % enable super- and sub-script notation

% Title Size & Weight
set(groot,'DefaultAxesTitleFontSizeMultiplier',1.4);
set(groot,'DefaultAxesTitleFontWeight','bold');
set(groot,'DefaultAxesLabelFontSizeMultiplier',1.2);

% Line Width, Color, and Style
if ~doFullColor
    colorSet = (0:.2:.6)'*ones(1,3);
    set(groot,'DefaultAxesColorOrder',colorSet);
end
set(groot,'DefaultAxesLineStyleOrder','-|--|:'); % Iterate over line styles after colors exhausted
set(groot,'DefaultLineLineWidth',1);

% Background Color
if doFullColor
    bgColor = .85*ones(1,3);
else
    bgColor = .93*ones(1,3);
end
set(groot,'DefaultAxesColor',bgColor);
set(groot,'DefaultFigureColor',[1 1 1]); % set white for outside the axes

% Grid
set(groot,'DefaultAxesXGrid','on','DefaultAxesYGrid','on');
set(groot,'DefaultAxesGridColor','w');
set(groot,'DefaultAxesMinorGridColor','w');
set(groot,'DefaultAxesGridAlpha',1);
set(groot,'DefaultAxesLineWidth',.25);

% Figure Size - set to 6 inches wide by 4 inches tall,
% use a .25 inch border all around
posVec = [0 0 6 4];
set(groot,'DefaultAxesUnits','inches');
set(groot,'DefaultAxesOuterPosition',posVec);
set(groot,'DefaultFigurePosition',posVec);
set(groot,'DefaultFigureUnits','inches');
%set(groot,'DefaultFigureOuterPosition',posVec);

