function create_counterfactual_scenarios
%% Import baseline scenario

global c

weights = c.filter_weights;

load([c.model_folder, 'alg_inputs.mat']);


%% Decompose shocks if needed

if c.sh == 0 % actual
    z = baseline.z;
elseif c.sh == 1 % no sectoral shocks
    [z, ~] = decompose_shocks(baseline.z, weights);
elseif c.sh == 2 % no sectoral and residual shocks
    [~, z] = decompose_shocks(baseline.z, weights);
end


%% Trade costs
[n_countries, ~, i_services, n_years] = size(baseline.kappa);

if c.tc == 0 % actual
    kappa = baseline.kappa;
elseif c.tc == 1 % 1972
    kappa = repmat(baseline.kappa(:, :, :, 1), [1 1 1 n_years]);
elseif c.tc == 2 % free trade
    kappa = ones(size(baseline.kappa));
    kappa(:, :, i_services, :) = ...
        repmat(eye(n_countries), [1, 1, 1, n_years]);    
end


%% Save
baseline.z = z;
baseline.kappa = kappa;

save([c.model_folder, 'alg_inputs.mat'], 'baseline')

end