function [ ret ] = variable_T_by_delta_T( time_number, delta_T, var_omega, ...
                                          physic_grid, front_coordinate, ...
                                          angle_of_T_gradient_solid, angle_of_T_gradient_liquid )

one_T = line_by_x_and_angle( front_coordinate, angle_of_T_gradient_liquid, 1 );

delta_omega = abs(atan(front_coordinate/(one_T - delta_T/2)) - atan(front_coordinate/(one_T + delta_T/2)));
normalization_coefficient = delta_omega/angle_of_T_gradient_liquid;

temp = variable_T_distribution( time_number, normalization_coefficient, var_omega, ...
                                physic_grid, front_coordinate, ...
                                angle_of_T_gradient_solid, angle_of_T_gradient_liquid );
ret = [temp(1), temp(end)];

end

