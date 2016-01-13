function create_counterfactual_scenarios
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
[z_no_sector, z_no_sector_residual] = decompose_shocks(z, weights);



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


%% Save different scenarios
scenario = 'No sectoral productivity shocks';
scenario_id = 1;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector, kappa);

scenario = 'No sectoral and residual productivity shocks';
scenario_id = 2;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector_residual, kappa);

scenario = '1972 trade costs';
scenario_id = 3;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z, kappa_initial);

scenario = 'No sectoral productivity shocks and 1972 trade costs';
scenario_id = 4;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector, kappa_initial);

scenario = 'No sectoral and residual productivity shocks and 1972 trade costs';
scenario_id = 5;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector_residual, kappa_initial);
               
scenario = 'Free trade in agriculture and manufacturing';
scenario_id = 6;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z, kappa_free_trade);

scenario = 'No sectoral productivity shocks and free trade in agriculture and manufacturing';
scenario_id = 7;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector, kappa_free_trade);               

scenario = 'No sectoral and residual productivity shocks and free trade in agriculture and manufacturing';
scenario_id = 8;
save_counterfactual(out_folder, theta, scenario, scenario_id,...
                   z_no_sector_residual, kappa_free_trade);               
end