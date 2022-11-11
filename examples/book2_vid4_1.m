% Code demonstrated in Video 4.1 of the second texbook

% Coordinates for the Lincoln Memorial
ref_lla = [38.88939, -77.05006, 0];

[x,y,z] = utils.lla2ecef(ref_lla(1), ref_lla(2), ref_lla(3),'deg','m');
ref_ecef = [x, y, z];

% Coordinates for the Capitol Building
test_lla = [38.89017, -77.00909, 0];

[x,y,z] = utils.lla2ecef(test_lla(1), test_lla(2), test_lla(3),'deg','m');
test_ecef = [x, y, z];

% Compute distance via ECEF
dist_km = norm(test_ecef-ref_ecef)/1e3;
fprintf('Distance from Lincoln Memorial to Capitol Building is %.2f km.\n',dist_km);

% Compare to Google Maps approx.
dist_gmap_mi = 2.3;
dist_gmap_km = dist_gmap_mi * unitsratio('m','mi')/1e3;
fprintf('Distance estimated via Google Maps is %.2f km.\n',dist_gmap_km);

% Compute distance via ENU
[e,n,u] = utils.lla2enu(test_lla(1), test_lla(2), test_lla(3),...
                        ref_lla(1), ref_lla(2), ref_lla(3), 'deg', 'm');
test_enu = [e, n, u];
dist_enu_km = norm(test_enu)/1e3;
fprintf('Distance computed via ENU coordinates is %.2f km.\n',dist_enu_km);
fprintf('Error between ECEF/ENU distance calculations is: %.4f m.\n', abs(dist_km - dist_enu_km)*1e3);

% Compute distance via AER
[a,e,r] = utils.enu2aer(test_enu(1), test_enu(2), test_enu(3), 'm');
test_aer = [a, e, r];

fprintf('Distance computed via AER coordinates is %.2f km.\n',test_aer(3)/1e3);
fprintf('Error between ECEF/AER distance calculations is: %.4f m.\n', abs(test_aer(3) - dist_km*1e3));

fprintf('Bearing from Lincoln Memorial to U.S. Capitol Building is %.2f deg.\n',test_aer(1));
fprintf('Elevation angle of U.S. Capitol as seen by Lincoln Memorial is %.4f deg.\n', test_aer(2));

% Visualize LOS
[Z,R] = readgeoraster('nasadem_washingtondc.tif','OutputType','double');
los2(Z,R,test_lla(1), test_lla(2), ref_lla(1), ref_lla(2), 0);
grid on;
title('Line of Sight from Lincoln Memorial (left) to U.S. Capitol Building (right)')

los2(Z,R,test_lla(1), test_lla(2), ref_lla(1), ref_lla(2), 10);
grid on;
title('Line of Sight from 10 ft above ground at Lincoln Memorial (left) to U.S. Capitol Building (right)')
