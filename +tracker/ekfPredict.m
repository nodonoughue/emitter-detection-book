function [x_pred, P_pred] = ekfPredict(x_est, P_est, Q, f_fun, g_fun)

x_pred = f_fun(x_est);

F = g_fun(x_est);
P_pred = F * P_est * F' + Q;