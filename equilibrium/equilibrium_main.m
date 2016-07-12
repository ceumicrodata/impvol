function equilibrium_main

global c

theta = c.theta;

include_counterfactuals = c.include_counterfactuals;

in_folder = c.data_folder_algorithm_input;
out_folder = c.results_folder;
in_file = [in_folder, 'data_theta_', num2str(theta), '.mat'];
load(in_file)

equilibrium = struct('scenario', {},...
                     'P_nt', {}, 'P_njt', {},...
                     'w_nt', {}, 'w_njt', {},...
                     'L_nt', {}, 'L_njt', {},...
                     'd', {});

%%
% Calculate baseline equilibrium
% equilibrium(1) = equilibrium_algorithm(baseline);

% out_file = [out_folder, 'equilibrium_baseline_theta_', num2str(theta), '.mat'];
% save(out_file, 'equilibrium')


%%
if include_counterfactuals
    
    % Calculate counterfactual equilibria
    n_counterfactuals = length(counterfactual);
    
    for cf = 3:n_counterfactuals
        equilibrium_input = baseline;
        
        % Overwrite the inputs corresponding to the actual counterfactual
        equilibrium_input.scenario = counterfactual(cf).scenario;
        equilibrium_input.z = counterfactual(cf).z;
        equilibrium_input.kappa = counterfactual(cf).kappa;
        
        
        counterfactual_equilibrium = equilibrium_algorithm(equilibrium_input);
        
        
        equilibrium(cf + 1) = counterfactual_equilibrium;
        
        out_file = [out_folder, sprintf('equilibrium_%d_theta_', cf),...
            num2str(theta), '.mat'];
        save(out_file, 'counterfactual_equilibrium')
    end

end % if include_counterfactuals

out_file = [out_folder, 'equilibrium_theta_', num2str(theta), '.mat'];
save(out_file, 'equilibrium')
    
