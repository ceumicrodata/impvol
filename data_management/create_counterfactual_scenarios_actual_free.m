function create_counterfactual_scenarios_actual_free
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
% [Z_no_sector, Z_no_sector_residual] = decomposeshocks(Z, weights);



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


%% Interpolated trade costs

steps = 10;
gamma = [1:-(1 / steps):0];

for i = 1:(steps + 1)      
    scenario = sprintf('%.1f * actual costs   +   %.1f * free trade ', gamma(i), 1 - gamma(i));
    scenarioId = i;
    
    kappa_actual_step = gamma(i) * kappa + (1 - gamma(i)) * kappa_free_trade;
    
    savecounterfactual(outFolder, theta, scenario, scenarioId,...
                       Z, kappa_actual_step);
                 
end
               
                
end