function x = angle_bisector(x_sensor,psi)
% x = angle_bisector(x_sensor,psi)
%
% Compute the center via intersection of angle bisectors for  
% 3 or more LOBs given by sensor positions xi and angle 
% of arrival measurements psi.
%
% If more than 3 measurements are provided, this method repeats on each set
% of 3 measurements, and then takes the average of the solutions.
%
% INPUTS:
%
%   x_sensor    N x 2 matrix of sensor positions
%   psi         N x 1 vector of AOA measurements (radians)
%
% OUTPUTS:
%
%   x           1 x 2 vector of estimated source position
%
% Nicholas O'Donoughue
% 1 July 2019

%% Set up all possible permutations of 3 elements
N = numel(psi);
assert(N>2,'At least 3 DF measurements required for angle bisector.');

if N <= 3
    A = 1:N;
elseif N <= 15
    A = nchoosek(1:N,3);
else
    % If there are more than ~15 rows, nchoosek returns an extremely large
    % matrix. Beyond that point, it may be better to include each
    % measurement only once, rathern than computing all possible sets of 3
    N = N - mod(N,3); % Ignore the last 1 or 2 sensors
    
    % Generate a matrix A with the different indices
    A = reshape(1:N,[],3);
end

%% Loop over sets
N_sets = size(A,1);
x_est = zeros(2,N_sets);
for idx_set = 1:N_sets
    thisX = x_sensor(:,A(idx_set,:));
    thisPsi = psi(A(idx_set,:));
    
    % Find vertices
    v1 = utils.find_intersect(thisX(:,1),thisPsi(1),thisX(:,2),thisPsi(2));
    v2 = utils.find_intersect(thisX(:,2),thisPsi(2),thisX(:,3),thisPsi(3));
    v3 = utils.find_intersect(thisX(:,3),thisPsi(3),thisX(:,1),thisPsi(1));

    % Find angle bisectors
    VV = [v1;v2;v3];
    th_fwd = atan2(circshift(VV(:,2),-1,1)-VV(:,2),...
                   circshift(VV(:,1),-1,1)-VV(:,1));
    th_back = atan2(circshift(VV(:,2),1,1)-VV(:,2),...
                    circshift(VV(:,1),1,1)-VV(:,1));
    th_bisector = .5*(th_fwd+th_back);
    flipped = th_bisector < th_fwd;
    th_bisector(flipped) = th_bisector(flipped)+pi;

    % Find the intersection
    this_cntr = utils.find_intersect(v1,th_bisector(1),v2,th_bisector(2));
    x_est(:,idx_set) = this_cntr;
end

x = mean(x_est,2);