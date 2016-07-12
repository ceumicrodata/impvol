function analysis

global c

theta = c.theta;
weights = c.filter_weights;
in_folder = c.results_folder;
out_folder = c.results_folder;

in_file = [in_folder, 'equilibrium_theta_', num2str(theta), '.mat'];
load(in_file)

nScenarios = length(equilibrium);
n_countries = size(equilibrium(1).L_nt, 1);

variances = zeros(n_countries, nScenarios);

% load data
for scenarioID = 1:nScenarios
    valueAdded = equilibrium(scenarioID).L_nt .* equilibrium(scenarioID).w_nt;
    
    P_nt = equilibrium(scenarioID).P_nt;
    
    realGDP = valueAdded ./ P_nt;
    [~, cycle] = detrend_series(log(realGDP), weights);
    variances(:, scenarioID) = var(cycle, 0, 2);
end

out_file = [out_folder, 'variances_', num2str(theta), '.mat'];
save(out_file, 'variances')


in_file = [in_folder, 'data_rgdp_and_volatility.mat'];
load(in_file)

mean(data_volatility_total)
mean(variances)

% h = figure();
% bar([data_volatility_total, variances])
% title('Log Real GDP Volatility')
% xlim([0.5 25.5])
% xlabel('Country')
% legend('data', 'model', 'Location', 'NE')
% print(h, [out_folder, 'volatilities'], '-dpng', '-r0')