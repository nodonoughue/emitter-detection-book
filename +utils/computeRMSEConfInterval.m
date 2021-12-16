function confInterval = computeRMSEConfInterval(gamma)
% confInterval = computeRMSEConfInterval(gamma)
%
% Determines the confidence interval for a given scale factor gamma, which
% is defined as the percentage of a standard normal distribution that falls
% within the bounds -gamma to gamma.
%
% Computed simply with:
%   confInterval = normcdf(gamma) - normcdf(-gamma)
%
% Inputs:
%
%   gamma           Scale factor
%
% Outputs:
%
%   confInterval    Confidence interval on scale (0-1)
%
% Nicholas O'Donoughue
% 1 July 2019

if license('test','Statistics_Toolbox')
    confInterval = normcdf(gamma) - normcdf(-gamma);
else
    confInterval = utils.normcdf(gamma) - utils.normcdf(-gamma);
end