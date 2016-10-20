function run_model(model)
% RUN_MODEL Runs model with given specification.
% Steps:
% 1. Initialize global parameters
% 2. Calibrate shocks from data
% 3. Compute counterfactual trade costs and/or shocks (if needed)
% 4. Compute equilibrium

addpath(genpath(pwd))

%% Step 1.
global c
c = init_globals(model);


%% Step 2.
calibrate_shocks;


%% Step 3.
create_counterfactual_scenarios


%% Step 4.
equilibrium_main



%% Step 5.
compute_volatilities