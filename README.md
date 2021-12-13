# MATLAB Companion to Emitter Detection and Geolocation for Electronic Warfare

![Textbook cover](https://us.artechhouse.com/Assets/ProductImages/5645_t.jpg "Emitter Detection and Geolocation for Electronic Warfare")

This repository contains the MATLAB code, released under the MIT License, that generates all of the figures and implements all of the algorithms and many of the performance calculations within the text *Emitter Detection and Geolocation for Electronic Warfare,* by Nicholas A. O'Donoughue, Artech House, 2019.

The textbook can be purchased from **[Artech House](https://us.artechhouse.com/Emitter-Detection-and-Geolocation-for-Electronic-Warfare-P2052.aspx)**.

### Update 20211213 -- Additional Features Coming
This repository now contains code for several new features, that are to be discussed in an upcoming textbook.  When it is published, I will include details on how to purchase the book.  Until then, consider these new features to be unfinished, but please feel free to send me bug reports if you have any.

A few of the new features include:
+ Scripts for handling 3D coordinate systems (including ENU, ECEF, and Geodetic Lat/Lon)
+ The ability to specify arbitrary test/reference sensor pairs for TDOA and FDOA (as opposed to setting the final sensor to be a common reference), along with code to modify covariance matrices to translate from sensor-level errors (TOA/FOA on each sensor) to measurement-level (covariance of TDOA and FDOA measurements between sensor pairs)
+ Generic solvers for constrained cases


## Installation

Simply open the folder in MATLAB or add the root folder containing this repository to the MATLAB search path.  Clone this repository, or download the latest release from the Mathworks File Exchange.

[![View Emitter Detection and Geolocation for Electronic Warfare on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74044-emitter-detection-and-geolocation-for-electronic-warfare)

## Figures

To generate the figures in the book, execute the script **runAll.m**

*Warning: Some of the figures are generated through a Monte Carlo simulation that can take several hours to complete.  For this reason, the variable **force_recalc** is set to False, which will skip any figures that were generated with a lengthy simulation.  To generate those figures, set force_recalc=True.*

Alternatively, one can call one of the chapter specific scripts in the folder make_figures, such as:

    force_recalc=True;
	make_figures/ch1_drawFigures;

The output will be stored in the **figures/** folder.

## Examples
The **examples/** folder contains the code to execute each of the examples in the textbook.

## Homework
The **hw/** folder contains data sets used for two homework problems in Chapter 8.

## Live Scripts
The **live_scripts/** folder contains MATLAB Live script versions for the examples in **examples/**.  It is currently a work in progress, and more will be added over the coming weeks.  [MATLAB Live scripts](https://www.mathworks.com/help/matlab/matlab_prog/what-is-a-live-script-or-function.html) blend text and code to enable interactive learning.

## Utilities
A number of utilities are provided in this repository, under the following namespaces:

+ **aoa/** Code to execute angle-of-arrival estimation, as discussed in Chapter 7
+ **array/** Code to execute array-based angle-of-arrival estimation, as discussed in Chapter 8
+ **atm/** Code to model atmospheric loss, as discussed in Appendix Carlo
+ **detector/** Code to model detection performance, as discussed in Chapter 3-4
+ **fdoa/** Code to execute Frequency Difference of Arrival (FDOA) geolocation processing, as discussed in Chapter 12.
+ **hybrid/** Code to execute hybrid geolocation processing, as discussed in Chapter 13.
+ **noise/** Code to model noise power, as discussed in Appendix D.
+ **prop/** Code to model propagation losses, as discussed in Appendix B.
+ **tdoa/** Code to execute Time Difference of Arrival (TDOA) geolocation processing, as discussed in Chapter 11.
+ **triang/** Code to model triangulation from multiple AOA measurements, as discussed in Chapter 10.
+ **utils/** Generic utilities, including numerical solvers used in geolocation algorithms.

## Feedback
Please submit any suggestions, bugs, or comments to nicholas [dot] odonoughue [at] ieee [dot] org.
