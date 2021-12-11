function [x_est,A,x_grid] = bfSoln(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,epsilon,ref_idx,pdftype)
% [x_est,A,x_grid] = bfSoln(x_fdoa,v_fdoa,rho_dot,C,x_ctr,search_size,epsilon,ref_idx,pdfType)
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
%   x_fdoa      Sensor positions [m]
%   v_fdoa      Sensor velocities [m/s]
%   rho_dot     Measurement vector [Hz]
%   C           Measurement error covariance matrix
%   x_ctr       Center of search grid [m]
%   search_size 2-D vector of search grid sizes [m]
%   epsilon     Desired resolution of search grid [m]
%   ref_idx     Scalar index of reference sensor, or nDim x nPair
%               matrix of sensor pairings
%   pdftype     String indicating the type of distribution to use.
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

if nargin < 8 || ~exist('ref_idx','var')
    ref_idx = [];
end

if nargin < 9 || isempty(pdftype)
    pdftype = [];
end

% Resample covariance matrix
n_sensor = size(x_fdoa, 2);
[test_idx_vec, ref_idx_vec] = utils.parseReferenceSensor(ref_idx, n_sensor);
C_tilde = utils.resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Generate the PDF
pdfs = utils.makePDFs(@(x) fdoa.measurement(x_fdoa,v_fdoa,x,ref_idx), ...
                     rho_dot,pdftype,C_tilde);

% Call the util function
[x_est,A,x_grid] = utils.bestfix(pdfs,x_ctr,search_size,epsilon);