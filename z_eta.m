function [ ret ] = z_eta( index, physic_grid, additional_border_nodes, eta_step )
% index - index of node for geting derivative
% physic_grid - real physical grid
% additional_border_nodes - pseudo pre-first and after-last elements for
%                           calculate derivarive with second approximation order
% eta_step - step by mathematic grid

calculation_grid = [additional_border_nodes(1), ...
                    physic_grid, ...
                    additional_border_nodes(2)];

if is_integer(index)
    % because we prepend virtual node
    real_index = index +1;
    ret = ( calculation_grid(real_index + 1) - ...
            calculation_grid(real_index - 1) )/(2*eta_step);
else
    % because we prepend virtual node
    real_index = floor(index) + 1;
    ret = ( calculation_grid(real_index + 1) - ...
            calculation_grid(real_index) )/eta_step;
end

end

