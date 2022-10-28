% Regenerate all figures from the textbook.
%
% WARNING: This will take a long time to run, some of the figures
%          rely on Monte Carlo simulations with a large number of
%          trials.
%
%          Set force_recalc=false to skip the figures that take
%          a long time.
%
% 25 July 2019
% Nicholas O'Donoughue

%% Emitter Detection and Geolocation Textbook
addpath('make_figures');

% Set this to true to run all of the Monte Carlo
% examples, although this may take several days to run.
force_recalc=false; %#ok<NASGU> 

ch1_drawFigures;
ch2_drawFigures;
ch3_drawFigures;
ch4_drawFigures;
ch5_drawFigures;
ch6_drawFigures;
ch7_drawFigures;
ch8_drawFigures;
ch9_drawFigures;
ch10_drawFigures;
ch11_drawFigures;
ch12_drawFigures;
ch13_drawFigures;
appB_drawFigures;
appC_drawFigures;
appD_drawFigures;

%% Practical Geolocation Textbook

% Set this to true to run all of the Monte Carlo
% examples, although this may take several days to run.
force_recalc=true;

book2_ch1_drawFigures;
book2_ch2_drawFigures;
book2_ch3_drawFigures;
book2_ch4_drawFigures;
book2_ch5_drawFigures;
book2_ch6_drawFigures;
book2_ch7_drawFigures;
book2_ch8_drawFigures;