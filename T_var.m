function [ ret ] = T_var( index, T_array, T_borders )
% index - index of element
% T_array - array with temperature
% T_borders - prefirst and afterlast nodes for temperature

real_T_array = [T_borders(1), T_array, T_borders(2)];
if is_integer(index)
    % not used
else
    real_index = index + 1/2;
    ret = real_T_array(real_index);
end

end

