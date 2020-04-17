function result = is_feasible(xtx,r_bistatic,res,min_dist)





% Defaults
if nargin < 4 || isempty(min_dist)
    min_dist = 1e3;
end

if nargin < 3 || isempty(res)
    res = min_dist/5;
end

% Build test vector
grid_max_xy = r_bistatic*.6;
grid_max_z =  20e3;
xy_vec = -grid_max_xy:res:grid_max_xy;
z_vec = 0:res:grid_max_z;

[x,y,z] = ndgrid(xy_vec,vy_vec,z_vec);
x_grid = [x,y,z];
Rt = utils.rng(xtx,x_grid);
Rr = utils.rng([0 0 0],x_grid);

min_err = min(abs(Rt + Rr - r_bistatic(:)'),[],2);

result = min(min_err) <= min_dist;