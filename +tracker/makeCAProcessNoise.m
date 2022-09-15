function Q = makeCAProcessNoise(maxG, num_states, heading, pitch, T, type,prob_no_accel,prob_max_accel)
% Q = makeCAProcessNoise(maxG, heading, pitch, T, type)
%
% Define the 9x9 process noise matrix for a 3D problem  with pos/vel/accel
% state space, with constant acceleration kinematic models, using one of
% three assumptions.
%
% type='lateral-only'
%       Air vehicle will divert laterally (x/y) with up to  maxG G's of 
%       lateral divert.  No thrusting changes.
%
% type='mostly-lateral'
%       Same as lateral only, except that 10% of divert process noise is
%       vertical, and 90% is lateral
%
% type='3d-divert'
%       Same as lateral only case, except the divert is allowed to
%       be in any dimension (x/y/z).
%
% To convert the loading factor to lateral acceleration, we use the 
% following equation:
%  max_accel = 9.81 * tan(acos(1./maxG))
%  accel_var = (1/3) * max_accel^2 * (1 + 4*prob_max_accel - prob_no_accel)
%
% INPUTS:
%   maxG        Maximum load factor (in G's)
%   heading     Heading (degrees CCW from +x axis)
%   pitch       Degrees above or below horizontal of target trajectory
%   T           Update rate (seconds)
%   type        String indicating the type of process noise to generate
%   prob_no_accel Probability that a given update will not experience a
%                 coordinated maneuver
%   prob_max_accel Probability that a coordinated maneuver will be
%                  performed at maxG
%
% OUTPUTS:
%   Q           9 x 9 covariance matrix for process noise
%
% Nicholas O'Donoughue
% 10 Feb 2022

if nargin < 7 || isempty(prob_no_accel)
    prob_no_accel = .7;
end
if nargin < 8 || isempty(prob_max_accel)
    prob_max_accel = .2;
end

% convert load factor to lateral acceleration term for level flight
max_accel = 9.81*tan(acos(1./maxG)); % m/s^2
accel_var = (1/3) * max_accel.^2 * (1 + 4*prob_max_accel - prob_no_accel); % convert inst. accel var to process noise power

if numel(accel_var) > 1
    accel_var = shiftdim(accel_var,-2); % shift any multi-dimensional inputs to the 3rd dimension and after
end
if numel(heading) > 1
    heading = shiftdim(heading,-2);
end
if numel(pitch) > 1
    pitch = shifdim(pitch,-2);
end

switch lower(type)
    case 'lateral-only'
        % No vertical acceleration, use heading along for lateral divert

        % 9x9 Process Noise Matrix
        row_x = cat(2,cosd(heading).^2, -cosd(heading).*sind(heading), zeros(size(heading)));
        row_y = cat(2,-cosd(heading).*sind(heading), sind(heading).^2, zeros(size(heading)));
        row_z = zeros(size(row_x));
        
        qq_xyz = accel_var .* cat(1,row_x,row_y,row_z);
    case '3d-divert'
        % No thrusting, only divert (in 3D)
        % Good for maneuvering aircraft under constant velocity
        ur = cat(2,cosd(heading).*cosd(pitch), sind(heading)*cosd(pitch), sind(pitch)+zeros(size(heading)));
        if numel(size(ur)<=2)
            uc = ur(:);
        else
            uc = permute(ur,[2,1,3:numel(size(ur))]);
        end
        P = eye(3) - uc.*conj(ur);
        qq_xyz = accel_var .* P;
%        qq_xyz = accel_std_dev^2/3 * [cosd(heading)^2*cosd(pitch)^2, -cosd(heading)*sind(heading)*cosd(pitch)^2, cosd(heading)*sind(pitch)*cosd(pitch);
%                                -cosd(heading)*sind(heading)*cosd(pitch)^2, sind(heading)^2*cosd(pitch)^2, sind(heading)*cosd(pitch)*sind(pitch);
%                                cosd(heading)*cosd(pitch)*sind(pitch), sind(heading)*cosd(pitch)*sind(pitch), sind(pitch)^2];
    case 'mostly-lateral'
        % Allow 10% of divert in vertical (z) dimension, 90% in lateral
        % 9x9 Process Noise Matrix
        row_x = cat(2,cosd(heading).^2, -cosd(heading).*sind(heading), zeros(size(heading)));
        row_y = cat(2,-cosd(heading).*sind(heading), sind(heading).^2, zeros(size(heading)));
        row_z = zeros(size(row_x));
        
        qq_xyz = accel_var .* cat(1,row_x,row_y,row_z);

        warning('mostly lateral case not properly defined');
end

% Pos/Vel/Accel distribution
j = [.5*T^2;T;1]; 
j = j(1:num_states); % keep only the first num_states states
J = j*j';

% Distribute x/y/z process noise to position, accel, and vel states
Q = kron(J,qq_xyz);