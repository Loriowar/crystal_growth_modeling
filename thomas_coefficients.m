function [ a, b, c, F ] = thomas_coefficients( physic_grid_previous, physic_grid_current, ...
                                               physic_grid_previous_time_layer, ...
                                               T_array, T_borders, ...
                                               k_star_index, ...
                                               alpha_sl, beta, delta, ...
                                               eta_step, tau_step, ...
                                               lambda_sl, stefan_const )
% generate coefficients for equation a_iT(i-1) - c_iT(i) + b_iT(i+1) = -F_i
% input parameters:
%  grid_params - [K, k*] - array with grid size and border coordinate
%  alpha - tau*mu/(2*h_eta) (array for solid and liquid: [alpha_s, alpha_l])
%  delta - tau/(4*h_eta)
%  h_eta - step by event coordinate grid

K = length(physic_grid_previous);
k_star = k_star_index;

% initialize variables for output coefficients
a = [];
b = [];
c = [];
F = [];

all_parts_coef = 1:K; % in document 0..K-1
border_nodes = [k_star-1, k_star, k_star+1]; % depricated
% generate array without SL border indexes
coef_without_borders_solid = all_parts_coef(2:k_star-2);
coef_without_borders_liquid = all_parts_coef(k_star+2:end-1);

% prepare function for simple assing coefficients
C_coef_simp = @(k_index) C_coef( k_index, k_star, ...
                                 physic_grid_previous, physic_grid_current, ...
                                 eta_step );
C_S_coef_simp = @(k_index) C_S_coef( k_index, k_star, ...
                                     physic_grid_previous, physic_grid_current, ...
                                     eta_step );
C_L_coef_simp = @(k_index) C_L_coef( k_index, k_star, ...
                                     physic_grid_previous, physic_grid_current, ...
                                     eta_step );
z_tau_simp = @(index)  z_tau( index, ...
                              physic_grid_previous, physic_grid_current, ...
                              tau_step );
psi_next = @(k_index) psi_var( k_index, physic_grid_current, eta_step );
psi_prev = @(k_index) psi_var( k_index, physic_grid_previous_time_layer, eta_step );
T_simp = @(T_index) T_var( T_index, T_array, T_borders );

% generate coefficients
% coefficients for solid part
alpha = alpha_sl(1);
for i = coef_without_borders_solid
    a(i) = alpha*C_S_coef_simp(i) - delta*z_tau_simp(i);
    b(i) = alpha*C_S_coef_simp(i+1) + delta*z_tau_simp(i+1);
    c(i) = alpha*C_S_coef_simp(i) + delta*z_tau_simp(i) + ...
           alpha*C_S_coef_simp(i+1) - delta*z_tau_simp(i+1) + psi_next(i+1/2);
    F(i) = psi_prev(i+1/2)*T_simp(i+1/2) + ...
           alpha*C_S_coef_simp(i+1)*(T_simp(i+3/2) - T_simp(i+1/2)) - ...
           alpha*C_S_coef_simp(i)*(T_simp(i+1/2) - T_simp(i-1/2)) + ...
           delta*z_tau_simp(i+1)*(T_simp(i+3/2) + T_simp(i+1/2)) - ...
           delta*z_tau_simp(i)*(T_simp(i+1/2) + T_simp(i-1/2));
end

% coefficients for liquid part
alpha = alpha_sl(2);
for i = coef_without_borders_liquid
    a(i) = alpha*C_L_coef_simp(i) - delta*z_tau_simp(i);
    b(i) = alpha*C_L_coef_simp(i+1) + delta*z_tau_simp(i+1);
    c(i) = alpha*C_L_coef_simp(i) + delta*z_tau_simp(i) + ...
           alpha*C_L_coef_simp(i+1) - delta*z_tau_simp(i+1) + psi_next(i+1/2);
    F(i) = psi_prev(i+1/2)*T_simp(i+1/2) + ...
           alpha*C_L_coef_simp(i+1)*(T_simp(i+3/2) - T_simp(i+1/2)) - ...
           alpha*C_L_coef_simp(i)*(T_simp(i+1/2) - T_simp(i-1/2)) + ...
           delta*z_tau_simp(i+1)*(T_simp(i+3/2) + T_simp(i+1/2)) - ...
           delta*z_tau_simp(i)*(T_simp(i+1/2) + T_simp(i-1/2));
end

% last node in solid
alpha = alpha_sl(1);
i = k_star-1;
a(i) = alpha*C_S_coef_simp(i) - delta*z_tau_simp(i);
b(i) = 1/(2*eta_step)*(T_k_star_var('next') + T_k_star_var('cur'));
c(i) = psi_next(i+1/2) + 2*alpha*C_S_coef_simp(k_star) + ...
       alpha*C_S_coef_simp(i) + delta*z_tau_simp(i);
F(i) = psi_prev(i+1/2)*T_simp(i+1/2) + ...
       2*alpha*C_S_coef_simp(k_star)*T_k_star_var('next') + ...
       2*alpha*C_S_coef_simp(k_star)*( T_k_star_var('cur') - T_simp(i+1/2)) - ...
       alpha*C_S_coef_simp(i)*(T_simp(i+1/2) - T_simp(i-1/2)) - ...
       delta*z_tau_simp(i)*(T_simp(i-1/2) + T_simp(i+1/2));

% first node in liquid
alpha = alpha_sl(3);
i = k_star+1;
a(i) = -1/(2*eta_step)*(T_k_star_var('next') +  T_k_star_var('cur'));
b(i) = alpha*C_L_coef_simp(i) + delta*z_tau_simp(i);
c(i) = psi_next(i - 1/2) + alpha*C_L_coef_simp(i) + 2*alpha*C_L_coef_simp(k_star) - ...
       delta*z_tau_simp(i);
F(i) = psi_prev(i - 1/2)*T_simp(i - 1/2) + ...
       2*alpha*C_L_coef_simp(k_star)*T_k_star_var('next') + ...
       alpha*C_L_coef_simp(i)*(T_simp(i+1/2) - T_simp(i - 1/2)) - ...
       2*alpha*C_L_coef_simp(k_star)*(T_simp(i - 1/2) - T_k_star_var('cur')) + ...
       delta*z_tau_simp(i)*(T_simp(i+1/2) + T_simp(i - 1/2));
       
% border node
alpha = 0; % unused, for universally form of code
i = k_star;
a(i) = lambda_sl*C_S_coef_simp(i)*beta;
b(i) = C_L_coef_simp(i)*beta;
c(i) = -stefan_const;
F(i) = -beta*(C_L_coef_simp(i) + ...
        lambda_sl*C_S_coef_simp(i))*(T_k_star_var('next') +  T_k_star_var('cur')) + ...
        lambda_sl*C_S_coef_simp(i)*beta*T_simp(i - 1/2) + ...
        C_L_coef_simp(i)*beta*T_simp(i+1/2);

% remove empty elements (first and last)
a(1) = [];
b(1) = [];
c(1) = [];
F(1) = [];

a(end) = [];
b(end) = [];
c(end) = [];
F(end) = [];
end

