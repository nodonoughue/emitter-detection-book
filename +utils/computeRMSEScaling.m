function gamma = computeRMSEScaling(confInterval)
% gamma = computeRMSEScaling(confInterval)
%
% Computes the RMSE scaling factor for the specified confidence
% interval (between 0 and 1).  Defined as the integral limits (-gamma to
% gamma) that contain the desired percentage of random samples from a
% standard normal distrubtion (mean = 0, standard deviation = 1).
%
% It is computed simply with:
%   gamma = norminv(.5 + confInterval/2);
%
% and returns a value gamma such that
%   normcdf(gamma) - normcdf(-gamma) = confInterval
%
% Inputs:
%   confInterval        Confidence interval, between 0 and 1
%
% Outputs:
%   gamma               Scale factor to apply for RMSE scaling
%
% Nicholas O'Donoughue
% 1 July 2019

% Test for input validity
if confInterval <= 0 || confInterval >= 1
    error('Input out of bounds.  Confidence Interval must be between 0 and 1');
end

% Compute Gamma
if license('test','Statistics_Toolbox')
    gamma = norminv(.5 + confInterval/2);
else
    gamma = utils.norminv(.5 + confInterval/2);
end