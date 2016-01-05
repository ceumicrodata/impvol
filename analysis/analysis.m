function analysis

global c

theta = c.theta;
weights = c.filterWeights;
inFolder = c.resultsFolder;
outFolder = c.resultsFolder;

inFile = [inFolder, 'equilibrium_theta_', num2str(theta), '.mat'];
load(inFile)

nScenarios = length(equilibrium);
nCountries = size(equilibrium(1).L_nt, 1);

variances = zeros(nCountries, nScenarios);

% load data
for scenarioID = 1:nScenarios
    valueAdded = equilibrium(scenarioID).L_nt .* equilibrium(scenarioID).w_nt;
    
    P_nt = equilibrium(scenarioID).P_nt;
    
    realGDP = valueAdded ./ P_nt;
    [~, cycle] = detrendseries(log(realGDP), weights);
    variances(:, scenarioID) = var(cycle, 0, 2);
end

outFile = [outFolder, 'variances_', num2str(theta), '.mat'];
save(outFile, 'variances')

