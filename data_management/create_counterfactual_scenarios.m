function create_counterfactual_scenarios
%% Import baseline scenario

global c

theta = c.theta;
weights = c.filterWeights;

inFolder = c.dataFolderAlgorithmInput;
outFolder = c.dataFolderAlgorithmInput;

inFile = [inFolder, 'data_theta_', num2str(theta), '.mat'];
load(inFile);


%% Decompose shocks
Z = baseline.Z;
[Z_no_sector, Z_no_sector_residual] = decomposeshocks(Z, weights);



%% Trade costs
kappa = baseline.kappa;
[nCountries, ~, nSectors, nYears] = size(kappa);
iServices = nSectors;

% Fix trade costs on their initial level
kappaInitial = repmat(kappa(:, :, :, 1), [1 1 1 nYears]);

% Free trade kappa
kappa_free_trade = ones(size(kappa));
kappa_free_trade(:, :, iServices, :) = ...
    repmat(eye(nCountries), [1, 1, 1, nYears]);


%% Save different scenarios
scenario = 'No sectoral productivity shocks';
scenarioId = 1;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector, kappa);

scenario = 'No sectoral and residual productivity shocks';
scenarioId = 2;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector_residual, kappa);

scenario = '1972 trade costs';
scenarioId = 3;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z, kappaInitial);

scenario = 'No sectoral productivity shocks and 1972 trade costs';
scenarioId = 4;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector, kappaInitial);

scenario = 'No sectoral and residual productivity shocks and 1972 trade costs';
scenarioId = 5;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector_residual, kappaInitial);
               
scenario = 'Free trade in agriculture and manufacturing';
scenarioId = 6;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z, kappa_free_trade);

scenario = 'No sectoral productivity shocks and free trade in agriculture and manufacturing';
scenarioId = 7;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector, kappa_free_trade);               

scenario = 'No sectoral and residual productivity shocks and free trade in agriculture and manufacturing';
scenarioId = 8;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z_no_sector_residual, kappa_free_trade);               
end