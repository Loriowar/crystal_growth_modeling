function [ ret ] = z_tau( index, physic_grid_previous, physic_grid_current, tau_step )
% index - index of node for geting derivative
% physic_grid_previous - real physical grid on previous time step
% physic_grid_current - real physical grid on current time step
% tau_step - step by mathematic grid

% NOTE derivative by time use only with j+1/2 coefficient, i.e. exactly
% between two time layers. Than is why we dont get j_coefficieny in input
% arguments

if is_integer(index)
    ret = ( physic_grid_current(index) - ...
            physic_grid_previous(index) )/tau_step;
else
    % not used
end

end

