function [ A, B, C, F ] = alternative_thomas_coefficients_for_concentration( physic_grid_current, ...
                                                                             C_array, C_borders, ...
                                                                             k_star_index, ...
                                                                             eta_step, tau_step, ...
                                                                             ReinoldsSchmitt_const)

grid_size = length(physic_grid_current);
eps = 0.5*tau_step/(ReinoldsSchmitt_const*eta_step^2);
for i = (k_star_index+1):grid_size
    eps1(i) = eps*eta_step/(physic_grid_current(i) - physic_grid_current(i - 1));
end
% eps1 = eps*eta_step/diff(physic_grid_current((k_star_index+1):end));
pskfp = eta_step/(physic_grid_current(k_star_index + 1) - physic_grid_current(k_star_index));
for i = (k_star_index+2):(grid_size - 1)
    psin(i) = 2*eta_step/(physic_grid_current(i + 1) - physic_grid_current(i - 1));
end
psin(k_star_index+1) = eta_step/(physic_grid_current(k_star_index+1) - physic_grid_current(k_star_index));
psin(grid_size) = eta_step/(physic_grid_current(end) - physic_grid_current(end - 1));
for i = (k_star_index + 2):(grid_size - 1) % NOTE honestly speaking, here must be loop until 'grid_size'
    e1 = eps1(i);
    pp = psin(i);
    pm = psin(i - 1);
    
    A(i - k_star_index) = e1*pm;
    B(i - k_star_index) = e1*pp;
    C(i - k_star_index) = 1 + e1*(pp + pm);
    F(i - k_star_index) = e1*pm*C_array(i - k_star_index - 1) + ...
                          e1*pp*C_array(i - k_star_index + 1) + ...
                          (1 - e1*(pp + pm)*C_array(i - k_star_index));
end
e1 = eps1(k_star_index + 1);
pp = psin(k_star_index + 1);
pm = pskfp;

A(1) = e1*pm;
B(1) = e1*pp;
C(1) = 1 + e1*(pp + 2*pm);
F(1) = e1*pp*C_array(3) + (1 - e1*(pp + 2*pm))*C_array(2);
% A(1) = A(2);
% B(1) = B(2);
% C(1) = C(2);
% F(1) = F(2);
end

