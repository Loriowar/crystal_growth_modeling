function [ ret ] = C_L_coef( k_index, k_star, ...
                             physic_grid_previous, physic_grid_current, ...
                             eta_step )
% function for calculate C_coefficient in liquid
% k_index - index by coordinate eta
% j_index - index by time tau. All time it equal to j+1/2.
% physic_grid_previous - real physical grid on previous time step
% physic_grid_current - real physical grid on current time step
% eta_step - step by mathematic grid

if k_index == length(physic_grid_previous)
    % first approximation order
    cur_z_eta_val = (physic_grid_current(k_index) - physic_grid_current(k_index - 1))/eta_step;
    prev_z_eta_val = (physic_grid_previous(k_index) - physic_grid_previous(k_index - 1))/eta_step;
    ret = (cur_z_eta_val + prev_z_eta_val)/2;
    ret = 1/ret;
    return
end

if k_index == k_star
    % first approximation order
    cur_z_eta_val = (physic_grid_current(k_star + 1) - physic_grid_current(k_star))/eta_step;
    prev_z_eta_val = (physic_grid_previous(k_star + 1) - physic_grid_previous(k_star))/eta_step;
    ret = (cur_z_eta_val + prev_z_eta_val)/2;
    ret = 1/ret;
    return
end

ret = 1/z_eta_half_time( k_index, ...
                         physic_grid_previous, physic_grid_current, ...
                         [NaN NaN],  [NaN NaN], ...
                         eta_step );

end

