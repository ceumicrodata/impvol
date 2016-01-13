function equilibrium_main

global c

theta = c.theta;

includeCounterfactuals = c.includeCounterfactuals;

inFolder = c.dataFolderAlgorithmInput;
outFolder = c.resultsFolder;
inFile = [inFolder, 'data_theta_', num2str(theta), '.mat'];
load(inFile)

equilibrium = struct('scenario', {},...
                     'P_nt', {}, 'P_njt', {},...
                     'w_nt', {}, 'w_njt', {},...
                     'L_nt', {}, 'L_njt', {},...
                     'd', {});

% Calculate baseline equilibrium
equilibrium(1) = equilibrium_algorithm(baseline);

outFile = [outFolder, 'equilibrium_baseline_theta_', num2str(theta), '.mat'];
save(outFile, 'equilibrium')

if includeCounterfactuals
    
    % Calculate counterfactual equilibria
    nCounterfactuals = length(counterfactual);
    
    for cf = 1:nCounterfactuals
        equilibriumInput = baseline;
        
        % Overwrite the inputs corresponding to the actual counterfactual
        equilibriumInput.scenario = counterfactual(cf).scenario;
        equilibriumInput.Z = counterfactual(cf).Z;
        equilibriumInput.kappa = counterfactual(cf).kappa;
        
        
        counterfactual_equilibrium = equilibrium_algorithm(equilibriumInput);
        
        
        equilibrium(cf + 1) = counterfactual_equilibrium;
        
        outFile = [outFolder, sprintf('equilibrium_%d_theta_', cf),...
            num2str(theta), '.mat'];
        save(outFile, 'counterfactual_equilibrium')
    end

end % if includeCounterfactuals

outFile = [outFolder, 'equilibrium_theta_', num2str(theta), '.mat'];
save(outFile, 'equilibrium')
    
