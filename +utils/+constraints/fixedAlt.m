function [a, a_grad] = fixedAlt(alt, type)
%FIXEDALT Summary of this function goes here
%   Detailed explanation goes here

%% Parse Inputs
if nargin < 2 || isempty(type)
    type='ellipse';
end

%% Determine which constraint to return
switch lower(type)
    case 'flat'
        a = @(x) fixedAltConstraintFlat(x, alt);
        a_grad = @(x) fixedAltGradFlat(x);

    case 'sphere'
        a = @(x) fixedAltConstraintSphere(x, alt);
        a_grad = @(x) fixedAltGradSphere(x);

    case 'ellipse'
        a = @(x) fixedAltConstraintEllipse(x, alt);
        a_grad = @(x) fixedAltGradEllipse(x);

    otherwise
        error('Invalid case type when calling fixedAlt; must be one of {flat|sphere|ellipse}.');
end

%% Subordinate Functions
%  These functions implement the specific constraints and gradients, and
%  are defined in this manner because they are too complex to be directly
%  defined as an anonymous function.

    function [epsilon, x_valid] = fixedAltConstraintFlat(x, alt)
        % Flat Earth; altitude is equivalent to the z-dimension
        epsilon = x(3,:) - alt;
        
        % To make x valid, replace x(3,:) with alt
        x_valid = x;
        x_valid(3,:) = alt;
    end

    function epsilon_grad = fixedAltGradFlat(x)
        n = size(x,2);
        epsilon_grad = cat(1,zeros(2,n), ones(1,n));
    end

    function [epsilon, x_valid] = fixedAltConstraintSphere(x, alt)
        % Implement equation 5.5, and the scale term defined in 5.9

        radius_tgt_sq = sum(abs(x).^2,1); % x'x
        epsilon = radius_tgt_sq - (utils.constants.radiusEarth + alt).^2; % eq 5.5
        scale = (utils.constants.radiusEarth + alt)/ sqrt(radius_tgt_sq); % eq 5.9, modified

        x_valid = scale * x; % origin is the center of the spherical Earth,
                             % just scale x
    end

    function epsilon_grad = fixedAltGradSphere(x)
        % Implements the gradient of equation 5.5. with respect to x.
        epsilon_grad = 2*x;
    end

    function [epsilon, x_valid] = fixedAltConstraintEllipse(x, alt)
        % Convert ECEF to LLA
        [lat, lon, h] = utils.ecef2lla(x(1,:), x(2,:), x(3,:));

        % Compare altitude to desired
        epsilon = h - alt;

        % Find the nearest valid x -- replace computed alt with desired
        [xx, yy, zz] = utils.lla2ecef(lat, lon, alt);
        x_valid = [xx(:), yy(:), zz(:)]';
    end

    function epsilon_grad = fixedAltGradEllipse(x)
        % Implements equations 5.18 through 5.26

        % Load constants
        e1sq = utils.constants.first_ecc_sq;
        a = utils.constants.semimajor_axis_km * 1e3;
        
        % Compute geodetic latitude
        [lat, ~, ~] = utils.ecef2lla(x(1,:), x(2,:), x(3,:));
        lat_rad = lat*pi/180;

        % Break position into x/y/z components
        xx = x(1,:);
        yy = x(2,:);
        zz = x(3,:);

        % Pre-compute some repeated terms
        xy_len_sq = xx.^2+yy.^2;
        xy_len = sqrt(xy_len_sq);
        zz_sq = zz.^2;

        sin_lat = sin(lat_rad);
        cos_lat = cos(lat_rad);

        % Compute gradient of geodetic latitude, equations 5.24-5.26
        dlat_dx = -xx.*zz*(1-e1sq) ./ (xy_len .* (zz.^2 + (1-e1sq)^2 * xy_len_sq));
        dlat_dy = -yy.*zz*(1-e1sq) ./ (xy_len .* (zz.^2 + (1-e1sq)^2 * xy_len_sq));
        dlat_dz = (1-e1sq)*xy_len ./ (zz_sq + (1-e1sq)*xy_len_sq);

        % Compute gradient of effective radius, equations 5.21-5.23
        dR_dx = a*e1sq*sin_lat.*cos_lat.*dlat_dx./(1-e1sq*sin_lat.^2).^(1.5);
        dR_dy = a*e1sq*sin_lat.*cos_lat.*dlat_dy./(1-e1sq*sin_lat.^2).^(1.5);
        dR_dz = a*e1sq*sin_lat.*cos_lat.*dlat_dz./(1-e1sq*sin_lat.^2).^(1.5);

        % Compute gradient of constraint (epsilon), equations 5.18-5.20
        de_dx = (xx - yy.^2./xx)./(cos_lat.*xy_len)-dR_dx;
        de_dy = (xx + yy)./(cos_lat.*xy_len)-dR_dy;
        de_dz = -dR_dz;

        epsilon_grad = [de_dx; de_dy; de_dz];

    end
end