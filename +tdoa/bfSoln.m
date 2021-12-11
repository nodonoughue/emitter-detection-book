function [x_est,A,x_grid] = bfSoln(x_sensor,psi,C,x_ctr,search_size,epsilon,ref_idx,pdftype)
% [x_est,A,x_grid] = bfSoln(x_sensor,psi,C,x_ctr,search_size,epsilon,pdftype)
%
% Construct the BestFix estimate by systematically evaluating the PDF at 
% a series of coordinates, and returning the index of the maximum.  
% Optionally returns the full set of evaluated coordinates, as well.
%
% Assumes a multi-variate Gaussian distribution with covariance matrix C,
% and unbiased estimates at each sensor.  Note that the BestFix algorithm
% implicitly assumes each measurement is independent, so any cross-terms in
% the covariance matrix C are ignored.
%
% INPUTS:
%   x_sensor    Sensor positions [m]
%   psi         Measurement vectors [m]
%   C           Measurement error covariance matrix [m^2]
%   x_ctr       Center of search grid [m]
%   search_size 2-D vector of search grid sizes [m]
%   epsilon     Desired resolution of search grid [m]
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   pdfType     String indicating the type of distribution to use.
%               See +utils/makePDFs.m for options.
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Ref:
%  Eric Hodson, "Method and arrangement for probabilistic determination of 
%  a target location," U.S. Patent US5045860A, 1990, 
%  https://patents.google.com/patent/US5045860A
%
%
% Nicholas A. O'Donoughue
% 8 August 2020

if nargin < 7 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 8 || isempty(pdftype)
    pdftype = [];
end

% Resample covariance matrix
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Generate the PDF
pdfs = utils.makePDFs(@(x) tdoa.measurement(x_sensor,x,ref_idx), ...
                     psi,pdftype,C_tilde);

% Call the util function
[x_est,A,x_grid] = utils.bestfix(pdfs,x_ctr,search_size,epsilon);