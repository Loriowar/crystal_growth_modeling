function [ ret ] = z_eta_half_time( index, physic_grid_previous, physic_grid_current, ...
                                           additional_border_nodes_prev, additional_border_nodes_cur, ... 
                                           eta_step )
% index - index of node for geting derivative
% physic_grid_previous - real physical grid on previous time step
% physic_grid_current - real physical grid on current time step
% additional_border_nodes - pseudo pre-first and after-last elements for
%                           calculate derivarive with second approximation order
% eta_step - step by mathematic grid

previous_time_z_eta = z_eta( index, physic_grid_previous, additional_border_nodes_prev, eta_step );
current_time_z_eta = z_eta( index, physic_grid_current, additional_border_nodes_cur, eta_step );

ret = (previous_time_z_eta + current_time_z_eta)/2;

end

