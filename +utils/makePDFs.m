function pdfs = makePDFs(msmt_function,msmts,pdftype,varargin)
% Generate a joint PDF or set of unitary PDFs representing the measurements
% 'msmts', given the measurement function handle 'msmt_function',
% covariance matrix 'C' and pdftype
%
% INPUTS:
%   msmt_function       A single function handle that will accept an nDim x
%                       nSource array of candidate emitter positions and
%                       return an nMsmt x nSource array of measurements
%                       that those emitters are expected to have generated.
%   msmts               The received measurements
%   pdftype             The type of distribution to assume
%   ... additional arguments as specified by PDF type ...
%
%
% OUTPUTS:
%   pdfs                Cell array of function handles, each of which
%                       accepts an nDim x nSource array of candidate source 
%                       positions, and returns a 1 x nSource array of 
%                       probabilities.
%
%
% Supported pdftype options:
%   'MVN'               Mutli-variate Normal.  Additional inputs:
%                       C       nMsmt x nMsmt covariance matrix
%   'Normal'            c       nMsmt x 1 array of individual variances
%
% Nicholas O'Donoughue
% 6 August 2020

if nargin < 3 || isempty(pdftype)
    pdftype = 'MVN';
end

switch pdftype
    case 'MVN'
        C = varargin{1};
        pdfs = {@(x) utils.mvnpdf(msmt_function(x)',msmts(:)',C)'};
            
    case 'Normal'
        c = varargin{1};
        pdfs = {@(x) utils.mvnpdf(msmt_function(x)',msmts(:)',diag(c))'};
            
    otherwise
        error('Unrecognized PDF type setting: %s.\n',pdftype);
end