function [ ret ] = line_by_x_and_angle( k_intersect_coord, alpha, x_var )
% k_star_coodr - coordinate of intersect Ox
% aplpha - angle in radial
% x_var - value of y in corresponding point

k = tan(alpha);
b = -k*k_intersect_coord;
line_func = @(x) k*x + b;

ret = line_func(x_var);

end

