function [ ret ] = variable_T_distribution( time_number, amplitude_scale_coefficient, var_omega, ...
                                            physic_grid, front_coordinate, ...
                                            angle_of_T_gradient_solid, angle_of_T_gradient_liquid )


% amplitude_scale_coefficient = 0.0001;
% var_omega = 0.2;
ret = generate_initial_T( physic_grid, ...
                          front_coordinate, ...
                          angle_of_T_gradient_solid, ...
                          angle_of_T_gradient_liquid*...
                            (1 + amplitude_scale_coefficient*sin(2*pi*var_omega*time_number)) );
end

