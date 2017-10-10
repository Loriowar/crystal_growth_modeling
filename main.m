function [  ] = main(  )

% main function for modeling growth crystall process

% initialize base contsnt
solid_nodes_count = 12;
liquid_nodes_count = 24;
% k_star_index = solid_nodes_count + 1;
initial_k_star = 0.1;
q = 1.05;
real_time_step = 0.96; % sec % 0.096
characteristic_growth_rate = 1.3*10^(-6);
initial_growth_speed = characteristic_growth_rate;
characteristic_length = 0.1; % meter
front_coordinate = initial_k_star;
% is tau dimension or dimensionless? (dimension time = 0.096sec)
tau_step = 1.25*10^(-5); % 0.0005*0.01*0.25 % 1.25*10^(-6)
eta_step = 1/(solid_nodes_count + liquid_nodes_count);
T_borders = [NaN, NaN]; % potentialy unused
angle_of_T_gradient_solid = pi/3;
angle_of_T_gradient_liquid = pi/12;
lambda_s = 1.3*10^(-5); % need to know (calc from Peclet const)
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
cycle_breaker = 10; % amount steps for stop pteration into time layer
finish_front_position = 0.2;
% unidealization params
var_nu_freq = 20; % Hz
var_omega = 1/var_nu_freq;
% var_omega = 0.0625;%0.006; % front temperature fluctuation frequency (go by iteration counter)
delta_T = 0.072*10^(-3); % temperature fluctuation (dimensionless) % 0.22*10^(-4)
front_coordinate_of_start_unideality = 0.15;
%constants for concentration
initial_C = 5*10^(-3);
D_const = 10^(-4); %cm^2/sec
ReinoldsSchmitt_const = 10;
k_seg_const = 1.9;
conc_alpha = tau_step*D_const/(2*eta_step^2);

% amplitude_scale_coefficient = 0.0001;
% var_omega = 0.2;
% function binds for simplicity
% variable_T_front = @(time_number, physic_grid, front_coordinate) ...
%                      variable_T_distribution( time_number, amplitude_scale_coefficient, var_omega, ...
%                                               physic_grid, front_coordinate, ...
%                                               angle_of_T_gradient_solid, angle_of_T_gradient_liquid );

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
iteration_delta_x_prev = iteration_delta_x_cur;

% array with growth speed on each iteration
growth_rate = initial_growth_speed;
% array with growth steps
x_steps = iteration_delta_x_cur;
time_iteration_counter = 1;
T_storage = initial_T;
physical_grid_storage = initial_physic_grid;
time_number = 0;
previous_physic_grid = initial_physic_grid;
% border_T_storage = [];

%concentration initial values
current_C = initial_C*ones(1, length(initial_physic_grid) - k_star_index + 1);
C_storage = current_C;
border_C_storage = current_C(1);

% front_coordinate = front_coordinate + iteration_delta_x_cur;
tic
dimensionless_time = tau_step;
while front_coordinate < finish_front_position
    time_number = time_number + 1;
    T_borders(1) = line_by_x_and_angle( front_coordinate, angle_of_T_gradient_solid, 0 );
    T_borders(2) = line_by_x_and_angle( front_coordinate, angle_of_T_gradient_liquid, 1 );
%     border_T_storage = [border_T_storage T_borders(1)];
    time_layer_count = 0;
    dimensionless_time = dimensionless_time + tau_step;
%     tic
    while true
        time_layer_count = time_layer_count + 1;
        front_coordinate_in_layer = front_coordinate + iteration_delta_x_cur;
        new_physic_grid = generate_crystal_grid( front_coordinate_in_layer, solid_nodes_count, liquid_nodes_count, q );
%         current_T = generate_initial_T( new_physic_grid, front_coordinate_in_layer, angle_of_T_gradient_solid, angle_of_T_gradient_liquid ); 
        % DOUBT: maybe we must recalculate first and last temperatures?
%         new_T = current_T; % from algorithm al last page of notes

        % TODO neet to ask where we must set border temperature: every time
        % inside time laler or one time after before iteration within time
%         T_borders(1) = line_by_x_and_angle( front_coordinate_in_layer, angle_of_T_gradient_solid, 0 ) ;%+ sin(2*pi*dimensionless_time)*1e-1;
%         T_borders(2) = line_by_x_and_angle( front_coordinate_in_layer, angle_of_T_gradient_liquid, 1 );

        % iteration within time layer
        [a, b, c, F] = thomas_coefficients( curent_physic_grid, new_physic_grid, ...
                                            curent_physic_grid, ...
                                            current_T, T_borders, ...
                                            k_star_index, ...
                                            alpha, beta, delta, ...
                                            eta_step, tau_step, ...
                                            lambda_sl, stefan_const); 

        var_kappa = zeros(1,2);
        % TODO ask about sequenses of T1 and T2 for mu var
        var_mu = T_borders;
        thomas_result = prog(a, b, c, F, var_kappa, var_mu);
    %     thomas_result(k_star_index)
%         [current_T', thomas_result' ]
%         current_T = thomas_result;
%         current_T = generate_initial_T( new_physic_grid, front_coordinate_in_layer, angle_of_T_gradient_solid, angle_of_T_gradient_liquid ); 
        iteration_delta_x_prev = iteration_delta_x_cur;
        iteration_delta_x_cur = thomas_result(k_star_index);
        
        previous_physic_grid = new_physic_grid;

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
%     toc
    growth_rate = [growth_rate, iteration_delta_x_cur/tau_step];
    x_steps = [x_steps, iteration_delta_x_cur];
    front_coordinate = front_coordinate + iteration_delta_x_cur
%     current_T = generate_initial_T( new_physic_grid, front_coordinate_in_layer, angle_of_T_gradient_solid*(1+0.0001*sin(2*pi*0.2*time_number)), angle_of_T_gradient_liquid ); 
    if front_coordinate > front_coordinate_of_start_unideality
        current_T = variable_T_front(time_number, new_physic_grid, front_coordinate);
    else
        current_T = generate_initial_T( new_physic_grid, front_coordinate_in_layer, angle_of_T_gradient_solid, angle_of_T_gradient_liquid ); 
    end
    curent_physic_grid = generate_crystal_grid( front_coordinate, solid_nodes_count, liquid_nodes_count, q );
    time_iteration_counter = [time_iteration_counter, time_layer_count];
    T_storage = [T_storage; thomas_result];
    physical_grid_storage = [physical_grid_storage; new_physic_grid];
    % WARNING maybe it is wrong way for set first iteration on new time layer
%     iteration_delta_x_cur = 0;
%     time_layer_count
%     iteration_delta_x_cur
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

figure
pseudo_x = 1:length(growth_rate);
plot(pseudo_x*real_time_step, growth_rate)
grid
title('Growth rate')
xlabel('t, sec')
ylabel('\nu, mm/hour')

figure
plot(pseudo_x*real_time_step, x_steps)
grid
title('Delta x')
xlabel('t, sec')
ylabel('\Delta x, dimensionless')

figure
stairs(pseudo_x, time_iteration_counter)
grid
title('Amount iretation within time layer')
xlabel('step number')
ylabel('N')

figure
plot((1:length(border_C_storage))*real_time_step, border_C_storage)
grid
title('Border concentration')
xlabel('t, sec')
ylabel('C, dimensionless')

figure
plot(cumsum(x_steps) + initial_k_star, border_C_storage)
grid
title('Border concentration')
temp_label = xlabel('$$\varphi$$');
set(temp_label,'Interpreter','latex')
ylabel('C, dimensionless')

figure
hold
grid
size_C = size(C_storage);
for i = 2:size_C(1)
    if mod(i, 100) == 0
        plot(1:size_C(2), C_storage(i, :))
    % 	pause(0.05)
    end
end
title('Concentration')

figure
hold
grid
size_T = size(T_storage);
for i = 2:size_T(1)
    if mod(i, 100) == 0
        plot( get_half_node_values( physical_grid_storage(i, :)), T_storage(i, :))
    %     pause(0.05)
    end
end
title('Temperature')

hold off

end

