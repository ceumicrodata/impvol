function save_counterfactual(out_folder, theta, scenario, scenario_id, z, kappa)
                        
in_file = [out_folder, 'data_theta_', num2str(theta), '.mat'];
load(in_file);                        
                        
counterfactual(scenario_id).scenario = scenario;
counterfactual(scenario_id).z = z;
counterfactual(scenario_id).kappa = kappa;

out_file = [out_folder, 'data_theta_', num2str(theta), '.mat'];
save(out_file, 'counterfactual', '-append');
end