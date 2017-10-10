function [ ret ] = line_by_two_points( x_coord, y_coord, x_val )

line_func = @(x) (x - x_coord(1))./(x_coord(2) - x_coord(1)) .* (y_coord(2) - y_coord(1)) + y_coord(1);
ret = line_func(x_val);

end

