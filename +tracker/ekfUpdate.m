function [xk, Pk] = ekfUpdate(x_prev, P_prev, zeta, C, z, H)
% 

%% Evaluate the Measurement and Jacobian at x_prev
this_z = z(x_prev);
this_H = H(x_prev);

%% Compute the Innovation (or Residual)
yk = zeta - this_z;

%% Compute the Innovation Covariance
Sk = this_H * P_prev * this_H' + C;

%% Compute the Kalman Gain
Kk = P_prev*this_H'/Sk;

%% Update the Estimate
xk = x_prev + Kk*yk;
Pk = (eye(size(P_prev))-Kk*this_H)*P_prev;