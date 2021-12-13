function [x_est,A,x_grid] = bfSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,search_size,epsilon,tdoa_ref_idx,fdoa_ref_idx,pdftype)
% [x_est,A,x_grid] = bfSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,search_size,epsilon,tdoa_ref_idx,fdoa_ref_idx,pdftype)
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
%   x_aoa           AOA sensor positions [m]
%   x_tdoa          TDOA sensor positions [m]
%   x_fdoa          FDOA sensor positions [m]
%   v_fdoa          FDOA sensor velocities [m/s]
%   zeta            Combined measurement vector
%   C               Combined measurement error covariance matrix
%   x_ctr           Center of search grid [m]
%   search_size     2-D vector of search grid sizes [m]
%   epsilon         Desired resolution of search grid [m]
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%   pdftype         String indicating the type of distribution to use.
%                   See +utils/makePDFs.m for options.
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

if nargin < 10 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 11 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 12 || isempty(pdftype)
    pdftype = [];
end

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = utils.parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = utils.parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = cat(2,tdoa_test_idx_vec, n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,tdoa_ref_idx_vec, n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_aoa = C(1:n_aoa, 1:n_aoa);
C_tfdoa = C(n_aoa+1:end, n_aoa+1:end);
C_tilde = blkdiag(C_aoa, utils.resampleCovMtx(C_tfdoa, test_idx_vec, ref_idx_vec));

% Generate the PDF
pdfs = utils.makePDFs(@(x) hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x, tdoa_ref_idx, fdoa_ref_idx), ...
                     zeta,pdftype,C_tilde);

% Call the util function
[x_est,A,x_grid] = utils.bestfix(pdfs,x_ctr,search_size,epsilon);