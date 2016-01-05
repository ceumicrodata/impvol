function savecounterfactual(outFolder, theta, scenario, scenarioId, Z, kappa)
                        
inFile = [outFolder, 'data_theta_', num2str(theta), '.mat'];
load(inFile);                        
                        
counterfactual(scenarioId).scenario = scenario;
counterfactual(scenarioId).Z = Z;
counterfactual(scenarioId).kappa = kappa;

outFile = [outFolder, 'data_theta_', num2str(theta), '.mat'];
save(outFile, 'counterfactual', '-append');
end