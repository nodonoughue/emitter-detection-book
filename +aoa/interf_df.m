function psi_est= interf_df(x1,x2,d_lam)
% psi_est= interf_df(x1,x2,d_lam)
%
% Compute the estimated angle of arrival for an interferometer, given the
% complex signal at each of two receivers, and the distance between them.
%
% Inputs:
%
%   x1          Signal vector from antenna 1
%   x2          Signal vector from antenna 2
%   d_lam       Antenna spacing, divided by the signal wavelength
%
% Outputs:
%
%   psi_est     Estimated angle of arrival [radians]
%
% Nicholas O'Donoughue
% 1 July 2019

% The inner product of the two signals is a sufficient statistic for the
% phase between them, in the presence of a single signal and Gaussian noise
y = x1(:)'*x2(:);

% Use atan2 to solve for the complex phase angle
phi_est = atan2(imag(y),real(y));

% Convert from phase angle to angle of arrival
psi_est = asin(phi_est/(2*pi*d_lam));