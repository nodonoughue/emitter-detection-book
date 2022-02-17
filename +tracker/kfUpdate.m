function [xk, Pk] = kfUpdate(x_prev, P_prev, zeta, C, H)
% 

%% Evaluate the Measurement and Jacobian at x_prev
this_z = H*x_prev;

%% Compute the Innovation (or Residual)
yk = zeta - this_z;

%% Compute the Innovation Covariance
Sk = H * P_prev * H' + C;

%% Compute the Kalman Gain
Kk = P_prev*H'/Sk;

%% Update the Estimate
xk = x_prev + Kk*yk;
Pk = (eye(size(P_prev))-Kk*H)*P_prev;