global global_growth_rate;
global global_x_steps;
global global_time_iteration_counter;
global global_physical_grid_storage;
global global_T_storage;

size_T = size(global_T_storage);
for i = 2:size_T(1)
   plot( get_half_node_values( global_physical_grid_storage(i, :)), global_T_storage(i, :))
   pause(0.05)
end