function create_counterfactual_scenarios_double_kappa
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
% [nCountries, ~, nSectors, nYears] = size(kappa);
% iServices = nSectors;

% Fix trade costs on their initial level
% kappaInitial = repmat(kappa(:, :, :, 1), [1 1 1 nYears]);

% Free trade kappa
% kappa_free_trade = ones(size(kappa));
% kappa_free_trade(:, :, iServices, :) = ...
%     repmat(eye(nCountries), [1, 1, 1, nYears]);


double_kappa = min(2 * kappa, 1);

%% Save different scenarios
scenario = 'Double kappa';
scenarioId = 1;
savecounterfactual(outFolder, theta, scenario, scenarioId,...
                   Z, double_kappa);
             
end