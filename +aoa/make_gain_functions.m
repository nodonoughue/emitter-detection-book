function [g,g_dot] = make_gain_functions(type,d_lam,psi_0)
% [g,g_dot] = make_gain_functions(type,d_lam,psi_0)
%
% Generate function handles for the gain pattern (g) and gradient (g_dot),
% given the specified aperture type, aperture size, and mechanical steering
% angle.
%
% Inputs:
%
%   type        String indicating the type of aperture requested.  Supports
%               'omni', 'adcock', and 'rectangulare'
%   d_lam       Aperture length, in units of wavelength (d/lambda)
%   psi_0       Mechanical steering angle (in radians) of the array
%               [default = 0]
%
% Outputs:
%
%   g           Function handle to the antenna pattern g(psi), for psi in
%               radians
%   g_dot       Function handle to the gradient of the antenna pattern,
%               g_dot(psi), for psi in radians.
%
% Nicholas O'Donoughue
% 1 July 2019

% type = 'Adcock' or 'Rectangular'
% params
%   d_lam = baseline (in wavelengths)
%   psi_0 = central angle

switch lower(type)
    case 'omni'
        g = @(psi) 1;
        g_dot = @(psi) 0;
    case 'adcock'
        g = @(psi) 2*sin(pi*d_lam.*cos(psi-psi_0));
        g_dot = @(psi) -2*pi*d_lam.*sin(psi-psi_0).*cos(pi*d_lam.*cos(psi-psi_0));
    case 'rectangular'
        g = @(psi) abs(sinc((psi-psi_0).*d_lam/pi)); % sinc includes implicit pi
        g_dot = @(psi) utils.sinc_deriv((psi-psi_0).*d_lam)*d_lam;
    otherwise
        error('Aperture type not supported.');
end
