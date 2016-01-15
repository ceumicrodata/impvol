function create_counterfactual_scenarios_2
%% Import baseline scenario

global c

theta = c.theta;
weights = c.filter_weights;

in_folder = c.data_folder_algorithm_input;
out_folder = c.data_folder_algorithm_input;

in_file = [in_folder, 'data_theta_', num2str(theta), '.mat'];
load(in_file);


%% Decompose shocks
z = baseline.z;
% [z_no_sector, z_no_sector_residual] = decompose_shocks(z, weights);



%% Trade costs
kappa = baseline.kappa;
[n_countries, ~, n_sectors, n_years] = size(kappa);
i_services = n_sectors;

% Fix trade costs on their initial level
kappa_initial = repmat(kappa(:, :, :, 1), [1 1 1 n_years]);

% Free trade kappa
kappa_free_trade = ones(size(kappa));
kappa_free_trade(:, :, i_services, :) = ...
    repmat(eye(n_countries), [1, 1, 1, n_years]);


%% Interpolated trade costs

steps = 10;
gamma = [1:-(1 / steps):0];

for i = 1:(steps + 1)      
    scenario = sprintf('%.1f * 1972 costs   +   %.1f * free trade ', gamma(i), 1 - gamma(i));
    scenario_id = i;
    
    kappa_actual_step = gamma(i) * kappa_initial + (1 - gamma(i)) * kappa_free_trade;
    
    save_counterfactual(out_folder, theta, scenario, scenario_id,...
                       z, kappa_actual_step);
end
               
                
end