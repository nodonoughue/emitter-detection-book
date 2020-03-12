function [zenith_loss, zenith_loss_o, zenith_loss_w] = calcZenithLoss(freq,alt_start,zenith_angle)
% zenith_loss, zenith_loss_o, zenith_loss_w = calcZenithLoss(freq,alt_start,zenith_angle)
%
% Computes the cumulative loss from alt_start [m] to zenith (100 km
% altitude), for the given frequencies (freq) in Hz and angle from zenith
% zenith_angle, in degrees.
%
% Does not account for refraction of the signal as it travels through the
% atmosphere; assumes a straight line propagation at the given zenith
% angle.
%
% INPUTS:
%   freq :          Carrier frequency [Hz]
%   alt_start :     Starting altitude [m]
%   zenith_angle :  Angle between line of sight and 
%                   zenith (straight up) [deg]
%
% OUTPUTS:
%   zenith_loss :   Cumulative loss to the edge of the atmosphere [dB]
%   zenith_loss_o : Cumulative loss due to dry air [dB]
%   zenith_loss_w : Cumulative loss due to water vapor [dB]
%
% 11 March 2020
% Nicholas O'Donoughue

if nargin < 2 || isempty(alt_start)
    alt_start = 0;
end

if nargin < 3 || isempty(zenith_angle)
    zenith_angle = 0;
end

%% Determine which dimension altitude layers should be created in
%  This is done to smoothly handle N-dimensional inputs that may or may
%  not be the same size (e.g. frequency varied on one dimension, altitude
%  on another).
n_dims_in_use = max([numel(size(freq)),numel(size(alt_start)),numel(size(zenith_angle))]);
layer_dim = n_dims_in_use+1;

%% Make Altitude Layers
% From ITU-R P.676-11(12/2017), layers should be set at exponential
% intervals
numLayers = 922; % Used for ceiling of 100 km
layerDelta = .0001*exp(((1:numLayers)-1)/100); % Layer thicknesses [km], eq 21
layerDelta = reshape(layerDelta,[ones(1,n_dims_in_use),numel(layerDelta)]);
layerTop = cumsum(layerDelta); % [km]
layerBottom = layerTop - layerDelta; % [km]
layerMid = (layerTop+layerBottom)/2;

% Drop layers below alt_start
alt_start_km = alt_start/1e3;
layerMask = layerTop >= min(alt_start_km);
layerBottom = layerBottom(layerMask);
layerMid = layerMid(layerMask);
layerTop = layerTop(layerMask);

% Lookup standard atmosphere for each band
atmStruct = atm.standardAtmosphere(layerMid*1e3);

% Compute loss coefficient for each band
[ao,aw] = atm.gasLossCoeff(freq,atmStruct.P,atmStruct.e,atmStruct.T);

% Account for off-nadir paths and partial layers
el_angle_deg = 90 - zenith_angle;
layerDeltaEff = utils.computeSlantRange(max(layerBottom,alt_start_km), layerTop, el_angle_deg, true);
layerDeltaEff(layerTop <= alt_start_km) = 0; % Set all layers below alt_start_km to zero

% Zenith Loss by Layer (loss to pass through each layer)
zenith_loss_by_layer_oxygen = ao.*layerDeltaEff;
zenith_loss_by_layer_water = aw.*layerDeltaEff;

% Cumulative Zenith Loss
zenith_loss_o = sum(zenith_loss_by_layer_oxygen,layer_dim); % Loss from ground to the bottom of each layer
zenith_loss_w = sum(zenith_loss_by_layer_water,layer_dim); % Loss from ground to the bottom of each layer
zenith_loss = zenith_loss_o + zenith_loss_w;

