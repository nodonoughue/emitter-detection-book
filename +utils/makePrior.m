function [prior, fim_prior] = makePrior(type, varargin)
% Return a prior distribution, and its Fisher Information Matrix. Currently
% only supports Gaussian (multivariate normal) as a valid type, for which
% the usage is:
%
% [prior, fim_prior] = makePrior('gaussian', mu, sigma)
%
% where mu is the expectation and sigma is the covariance matrix.
%
% INPUTS:
%   type                String dictating type of statistical prior;
%                       currently only 'gaussian' is implemented
%   mu                  n_dim x 1 expectation                
%   sigma               n_dim x n_dim covariance matrix
%
% OUTPUTS:
%   prior               Function handle to statistical prior; accepts an
%                       n_dim x N matrix of input positions
%   fim_prior           Function handle to compute the Fisher Information
%                       Matrix of the prior, at a given point.  Accepts 
%                       n_dim x N inputs, and returns an n_dim x n_dim x N 
%                       datacube of FIM for each input vector.
%
% Nicholas O'Donoughue
% 24 November 2021

switch lower(type)
    case {'gaussian','mvn'}
        assert(nargin ==3,'Not enough inputs.');

        mu = varargin{1};
        sigma = varargin{2};
        prior = @(x) utils.mvnpdf(x', mu', sigma);

        warning('fim_prior not yet implemented...DERIVE!');
        fim_prior = sigma;
    otherwise
        error('Unrecognized type setting: %s.\n',type);
end