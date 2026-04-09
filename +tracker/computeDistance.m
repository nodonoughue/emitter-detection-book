function [dist, innov, S] = computeDistance(s_pred, zeta, msmt_model)
% computeDistance  Compute Mahalanobis distance squared between a predicted
%                  state and a measurement.
%
% [dist, innov, S] = computeDistance(s_pred, zeta, msmt_model)
%
% INPUTS
%   s_pred      Predicted State struct (from predictState)
%   zeta        Measurement vector (num_msmt x 1)
%   msmt_model  Measurement model struct from makeMeasurementModel
%
% OUTPUTS
%   dist    Mahalanobis distance squared: y' * S^{-1} * y  (chi2-distributed with num_msmt DOF)
%   innov   Innovation vector y = zeta - z_hat  (num_msmt x 1)
%   S       Innovation covariance (num_msmt x num_msmt)
%
% Nicholas O'Donoughue
% June 2025

z_hat = msmt_model.z_fun(s_pred);
H     = msmt_model.h_fun(s_pred);

if ndims(H) == 3
    H = H(:, :, 1)';
end

innov = zeta(:) - z_hat(:);
S    = H * s_pred.covar * H' + msmt_model.R;
dist = innov' / S * innov;
