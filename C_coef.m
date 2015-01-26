function [ ret ] = C_coef( k_index, k_star, ...
                           physic_grid_previous, physic_grid_current, ...
                           eta_step )
% function for calculate C_coefficient in all parts of area, ecxept borders
% k_index - index by coordinate eta
% j_index - index by time tau. All time it equal to j+1/2.
% physic_grid_previous - real physical grid on previous time step
% physic_grid_current - real physical grid on current time step
% eta_step - step by mathematic grid

% set NaN for border nodex. It must be processed in other functions

if max(k_index == [1, length(physic_grid_previous), k_star])
    ret = NaN;
    return;
end

ret = 1/z_eta_half_time( k_index, ...
                         physic_grid_previous, physic_grid_current, ...
                         [NaN NaN],  [NaN NaN], ...
                         eta_step );
end

