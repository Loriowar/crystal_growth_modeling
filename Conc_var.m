function [ ret ] = Conc_var( index, C_array, C_borders, offset )
% index - index of element
% C_array - array with concentration
% C_borders - prefirst and afterlast nodes for concentration

real_C_array = [C_borders(1), C_array, C_borders(2)];
if is_integer(index)
    % TODO
else
    real_index = index + 1/2 - offset;
    ret = real_C_array(real_index);
end

end

