function [  ] = main(  )

% main function for modelling growth crystal process

% initialize base constants
solid_nodes_count = 24;
liquid_nodes_count = 48;
initial_k_star = 0.1;
q = 1.05;
real_time_step = 0.096; % sec
characteristic_growth_rate = 1.3*10^(-6);
initial_growth_speed = characteristic_growth_rate;
characteristic_length = 0.1; % meter
front_coordinate = initial_k_star;
tau_step = 1.25*10^(-6);
eta_step = 1/(solid_nodes_count + liquid_nodes_count);
T_borders = [NaN, NaN];
angle_of_T_gradient_solid = pi/3;
angle_of_T_gradient_liquid = pi/12;
lambda_s = 1.3*10^(-5);
lambda_l = 2.5*10^(-5);
peclet_number = characteristic_length * characteristic_growth_rate/lambda_s;
alpha_s = tau_step*lambda_s/(2*eta_step^2);
alpha_l = tau_step*lambda_l/(2*eta_step^2);
alpha_first_in_liquid = tau_step/(2*peclet_number*eta_step^2);
alpha = [alpha_s, alpha_l, alpha_first_in_liquid];
beta = tau_step/eta_step;
delta = tau_step/(4*eta_step);
lambda_sl = lambda_s/lambda_l;
stefan_const = 0.1185;
epsilon = 1e-12; % variable for stop iteration into time layer
cycle_breaker = 10; % amount steps for stop iteration into time layer
finish_front_position = 0.6;
% unidealization params
var_nu_freq = 1/10; % Hz
var_omega = 1/var_nu_freq;
delta_T = 0.22*10^(-4); % temperature fluctuation (dimensionless)
front_coordinate_of_start_unideality = 0.4;
%constants for concentration
initial_C = 5*10^(-5);
D_const = 10^(-4); %cm^2/sec
ReinoldsSchmitt_const = 10;
k_seg_const = 0.1;
conc_alpha = tau_step*D_const/(2*eta_step^2);

variable_T_front = @(time_number, physic_grid, front_coordinate) ...
                     variable_T_by_delta_T( time_number*real_time_step, delta_T, var_omega, ...
                                            physic_grid, front_coordinate, ...
                                            angle_of_T_gradient_solid, angle_of_T_gradient_liquid );

% initialize grid
[initial_physic_grid, k_star_index] = generate_crystal_grid( front_coordinate, ...
                                                             solid_nodes_count, liquid_nodes_count, q );
% NOTE: up and low temperatures set into function body
initial_T = generate_initial_T( initial_physic_grid, initial_k_star, angle_of_T_gradient_solid, angle_of_T_gradient_liquid ); 

% define variables for iteration process
curent_physic_grid = initial_physic_grid;
current_T = initial_T;

% delta x within iteration by time
iteration_delta_x_cur = initial_growth_speed*tau_step;

% array with growth speed on each iteration
growth_rate = initial_growth_speed;
% array with growth steps
x_steps = iteration_delta_x_cur;
time_iteration_counter = 1;
T_storage = initial_T;
physical_grid_storage = initial_physic_grid;
time_number = 0;
previous_physic_grid = initial_physic_grid;
new_physic_grid = previous_physic_grid;

%concentration initial values
current_C = initial_C*ones(1, length(initial_physic_grid) - k_star_index);
C_storage = current_C;
border_C_storage = current_C(1);

tic
dimensionless_time = tau_step;
while front_coordinate < finish_front_position
    time_number = time_number + 1;
    if front_coordinate > front_coordinate_of_start_unideality
        T_borders = variable_T_front(time_number, new_physic_grid, front_coordinate);
    else
        T_borders(1) = line_by_x_and_angle( front_coordinate, angle_of_T_gradient_solid, 0 );
        T_borders(2) = line_by_x_and_angle( front_coordinate, angle_of_T_gradient_liquid, 1 );
    end
    time_layer_count = 0;
    dimensionless_time = dimensionless_time + tau_step;
    while true
        time_layer_count = time_layer_count + 1;
        front_coordinate_in_layer = front_coordinate + iteration_delta_x_cur;
        new_physic_grid = generate_crystal_grid( front_coordinate_in_layer, solid_nodes_count, liquid_nodes_count, q );

        % iteration within time layer
        [a, b, c, F] = thomas_coefficients( curent_physic_grid, new_physic_grid, ...
                                            curent_physic_grid, ...
                                            current_T, T_borders, ...
                                            k_star_index, ...
                                            alpha, beta, delta, ...
                                            eta_step, tau_step, ...
                                            lambda_sl, stefan_const); 

        var_kappa = zeros(1,2);
        var_mu = T_borders;
        thomas_result = prog(a, b, c, F, var_kappa, var_mu);

        iteration_delta_x_prev = iteration_delta_x_cur;
        iteration_delta_x_cur = thomas_result(k_star_index);

        % normal condition for leaving time layer
        diff_norma = abs(iteration_delta_x_cur - iteration_delta_x_prev)/iteration_delta_x_cur;
        if time_layer_count > 1 && ...
           diff_norma < epsilon;
            break;
        end
        % protect from endless cycle
        if time_layer_count > cycle_breaker
            break
        end
    end
    
    %concentration calculate
    growth_rate_for_C = iteration_delta_x_cur/tau_step;
    new_Conc = calc_Conc_distribution( curent_physic_grid, new_physic_grid, ...
                                       curent_physic_grid, ...
                                       current_C, [NaN, NaN], ...
                                       k_star_index, ...
                                       conc_alpha, delta, ...
                                       eta_step, tau_step, ...
                                       growth_rate_for_C, ReinoldsSchmitt_const, k_seg_const);
    C_storage = [C_storage; new_Conc];
    border_C_storage = [border_C_storage, new_Conc(1)];
    current_C = new_Conc;
    
    current_T = thomas_result;

    growth_rate = [growth_rate, iteration_delta_x_cur/tau_step];
    x_steps = [x_steps, iteration_delta_x_cur];
    front_coordinate = front_coordinate + iteration_delta_x_cur
    curent_physic_grid = generate_crystal_grid( front_coordinate, solid_nodes_count, liquid_nodes_count, q );
    time_iteration_counter = [time_iteration_counter, time_layer_count];
    T_storage = [T_storage; thomas_result];
    physical_grid_storage = [physical_grid_storage; new_physic_grid];
end
toc

dimensionless_time
dimension_time = dimensionless_time*real_time_step/tau_step

global global_growth_rate;
global_growth_rate = growth_rate;
global global_x_steps;
global_x_steps = x_steps;
global global_time_iteration_counter;
global_time_iteration_counter = time_iteration_counter;
global global_physical_grid_storage;
global_physical_grid_storage = physical_grid_storage;
global global_T_storage;
global_T_storage = T_storage;
global global_C_storage;
global_C_storage = C_storage;

end

