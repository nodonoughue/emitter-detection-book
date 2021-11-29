function xy_lob = drawLob(x_sensor, psi, x_source, scale)
% function xy_lob = drawLob(x_sensor, psi, len)
%
% Draw a line of bearing starting at the coordinate x_sensor. Can handle
% multiple LOB drawings at once.
%
% INPUTS:
%   x_sensor        2xM matrix of sensor positions for M AOA sensors
%   psi             Mx1 vector of of AOA measurements (radians)
%   x_source        Optional 2x1 source position; used to compute baseline
%                   LOB distance for each measurement (default=1)
%   scale           Optional scalar multiple, use scale > 1 to draw LOB
%                   past the target.
%
% OUTPUTS:
%   xy_lob          2 x 2 x M datacube of LOB endpoints
%                   -- x coordinates in x(1,:,m) for the m-th LOB
%                   -- y coordinates in x(2,:,m) for the m-th LOB
%
% Nicholas O'Donoughue
% 27 April 2021

% Check inputs
M1 = size(x_sensor, 2);
M2 = numel(psi);

if M1 ~= M2
    fprintf('The number of sensor positions and measurements must match.\n');
    M = min(M1, M2);
    x_sensor = x_sensor(:,1:M);
    psi = psi(1:M);
else
    M=M1;
end

if nargin < 3 || isempty(x_source)
    % No source position, use unit range
    range = 1;
else
    % Compute range for each lob
    range = utils.rng(x_sensor, x_source);
end

if nargin < 4 || isempty(scale)
    % Default to unit scale
    scale = 1;
end

% Draw the LOBs
x_end = cos(psi(:)).*range*scale; % end points, x coordinate (Mx1)
y_end = sin(psi(:)).*range*scale; % end points, y coordinate (Mx1)

xy_end = cat(1,reshape(x_end,1,1,M), reshape(y_end,1,1,M)); % 2x1xM
xy_start = zeros(2,1,M); % 2x1xM
xy_lob_centered = cat(2,xy_start,xy_end);
xy_lob = reshape(x_sensor,2,1,M) + xy_lob_centered;