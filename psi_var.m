function [ ret ] = psi_var( k_index, physic_grid, eta_step )
% k_index - index by coordinate eta. It all time is half-integer.
% j_index - all time is integer and it calculate through physic_grid
% physic_grid - real physical grid

ret = z_eta( k_index, physic_grid, [NaN NaN], eta_step );

end

