function psi_est = watson_watt_df(r,x,y)
% psi_est = watson_watt_df(r,x,y)
%
% Compute the angle of arrival given a centrally located reference signal,
% and a pair of Adcock antennas oriented orthogonally.
%
% Inputs:
%
%   r           Reference signal from antenna centrally located between
%               the two Adcock receivers.
%   x           Test signal from primary Adcock receiver; oriented in the
%               +x direction (0 degrees)
%   y           Test signal from secondary Adcock reciever; oriented in the
%               +y direction (90 degrees CCW)
%
% Outputs:
%
%   psi_est     Estimated angle of arrival [radians]
%
% Nicholas O'Donoughue
% 1 July 2019

% Remove reference signal from test data
xx = (r(:)'*x(:));
yy = (r(:)'*y(:));

% Results should be V*cos(th) and V*sin(th), use atan2 to solve for th
psi_est = atan2(real(yy),real(xx)); % output in radians