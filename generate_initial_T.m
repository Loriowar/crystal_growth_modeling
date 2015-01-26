function [ ret ] = generate_initial_T( physic_grid, phisic_k_star, angle_s, angle_l )
% solid_nodes - amount nodes in solid part
% liquid_nodes - amount nodes in liquid part (must be even)

if length(phisic_k_star) > 1
    border_index = phisic_k_star(1);
    phisic_k_star = phisic_k_star(2);
else
    diff_border = 1e-6;
    for i = 1:length(physic_grid)
        if (physic_grid(i) + diff_border >= phisic_k_star && ...
            physic_grid(i) - diff_border <= phisic_k_star)
            border_index = i;
            break;
        end
    end
end

% temperature in left extremity of innerval (or in bottom)
down_T = line_by_x_and_angle( phisic_k_star, angle_s, 0 );

% temperature in right extremity of innerval (or in top)
up_T = line_by_x_and_angle( phisic_k_star, angle_l, 1 );

% border_index = find(physic_grid == phisic_k_star);
for i = 1:border_index - 1
    ret(i) = line_by_two_points([0, phisic_k_star], [down_T, 0], (physic_grid(i + 1) +  physic_grid(i))/2);
end

for i = border_index:length(physic_grid) - 1
    ret(i) = line_by_two_points([phisic_k_star, 1], [0, up_T], (physic_grid(i + 1) +  physic_grid(i))/2);
end

end

