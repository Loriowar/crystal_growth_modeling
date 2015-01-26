function [ half_nodes_grid ] = get_half_node_values( physical_grid )
half_nodes_grid = [];
for i = 2:length(physical_grid)
    half_nodes_grid = [half_nodes_grid, (physical_grid(i) + physical_grid(i - 1))/2];
end

end

