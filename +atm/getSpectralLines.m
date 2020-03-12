function [f_ctr_o,f_ctr_w] = getSpectralLines(fmin,fmax)
% [f_ctr_o,f_ctr_w] = getSpectralLines(fmin,fmax)
%
% Returns the frequencies for the spectral absorption lines of oxygen and
% water (in Hz).
%
% Inputs:
%   fmin        - Optional.  Minimum frequency to consider [Hz]
%   fmax        - Optional.  Maximum frequency to consider [Hz]
%
% Outputs:
%   f_ctr_o     - Array of spectral lines for oxygen absorption [Hz]
%   f_ctr_w     - Array of spectral lines for water vapor absorption [Hz]
%
% Nicholas O'Donoughue
% 1 July 2019

f_ctr_o = [50.474214,50.987745,51.50336,52.021429,52.542418,53.066934,53.595775,54.130025,54.67118,55.221384,55.783815,56.264774,56.363399,56.968211,57.612486,58.323877,58.446588,59.164204,59.590983,60.306056,60.434778,61.150562,61.800158,62.41122,62.486253,62.997984,63.568526,64.127775,64.67891,65.224078,65.764779,66.302096,66.836834,67.369601,67.900868,68.431006,68.960312,118.750334,368.498246,424.76302,487.249273,715.392902,773.83949,834.145546]*1e9;
f_ctr_w = [22.23508,67.80396,119.99594,183.310087,321.22563,325.152888,336.227764,380.197353,390.134508,437.346667,439.150807,443.018343,448.001085,470.888999,474.689092,488.490108,503.568532,504.482692,547.67644,552.02096,556.935985,620.700807,645.766085,658.00528,752.033113,841.051732,859.965698,899.303175,902.611085,906.205957,916.171582,923.112692,970.315022,987.926764,1780]*1e9;

if nargin >= 1 && ~isempty(fmin)
    f_ctr_o = f_ctr_o(f_ctr_o>=fmin);
    f_ctr_w = f_ctr_w(f_ctr_w>=fmin);
end

if nargin >= 2 && ~isempty(fmax)
    f_ctr_o = f_ctr_o(f_ctr_o<=fmax);
    f_ctr_w = f_ctr_w(f_ctr_w<=fmax);
end