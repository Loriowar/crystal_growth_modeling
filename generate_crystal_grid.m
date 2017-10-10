function [ ret, k_star_index ] = generate_crystal_grid( k_star, solid_nodes_count, liquid_nodes_count, q )
% k_star - index of front (belongs to [0, 1] interval)
% solid_nodes - amount nodes in solid part
% liquid_nodes - amount nodes in liquid part (must be even)
% q - geometric progression parameter  (1.01; 1.05; 1.10; 1.12)
solid_interval = [0, k_star];
liquid_interval = [k_star, 1];

% generate solid grid (eventual)
solid_step = (solid_interval(2) - solid_interval(1))/solid_nodes_count;
ret = cumsum(ones(1, solid_nodes_count)*solid_step);
% add begining of interval into result grid
ret = [0, ret];
% index of boundary in mathematic grid
k_star_index = solid_nodes_count + 1;

% generate liquid grid
half_liquid_nodes = liquid_nodes_count/2;
center_liquid_interval = (liquid_interval(2) + liquid_interval(1))/2;
half_liquid_interval_size = liquid_interval(2) - center_liquid_interval;

[~, first_liquid_grid_steps] = generate_unevent_grid( half_liquid_nodes, q );
first_liquid_grid_steps = first_liquid_grid_steps*(half_liquid_interval_size/liquid_interval(2));
liquid_grid_first_part_intervals = cumsum([ret(end),first_liquid_grid_steps]);
% remove first node, i.e. it is repeating node
liquid_grid_first_part_intervals(1) = [];
ret = [ret, liquid_grid_first_part_intervals];

last_liquid_grid_steps = fliplr(first_liquid_grid_steps);
liquid_grid_last_part_intervals = cumsum([ret(end),last_liquid_grid_steps]);
% remove first node, i.e. it is repeating node
liquid_grid_last_part_intervals(1) = [];
ret = [ret, liquid_grid_last_part_intervals];
end

