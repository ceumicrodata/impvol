%% Main file to replicate the impvol project
%% Preliminaries
clear all
clc
format compact
addpath(genpath(pwd))

%tic

global c

%% Initialize program parameters
c = initialize_io_ubt;



%% Data management
% Import and process data:
% Compute parameters and exogenous variables of the model.
data_management



%% Create counterfactual scenarios
% Compute alternative parameters and exogenous variables for the counterfactual
% scenarios.
create_counterfactual_scenarios



%% Model equilibrium
equilibrium_main



%% Volatility analysis
analysis



%% Model fit analysis
% model_fit_analysis



%% Main file to replicate the impvol project
%% Preliminaries
clear all
% clc
format compact
addpath(genpath(pwd))

%tic

global c

%% Initialize program parameters
c = initialize_io;



%% Data management
% Import and process data:
% Compute parameters and exogenous variables of the model.
data_management



%% Create counterfactual scenarios
% Compute alternative parameters and exogenous variables for the counterfactual
% scenarios.
create_counterfactual_scenarios



%% Model equilibrium
equilibrium_main



%% Volatility analysis
analysis




%% Main file to replicate the impvol project
%% Preliminaries
clear all
% clc
format compact
addpath(genpath(pwd))

%tic

global c

%% Initialize program parameters
c = initialize_ubt;



%% Data management
% Import and process data:
% Compute parameters and exogenous variables of the model.
data_management



%% Create counterfactual scenarios
% Compute alternative parameters and exogenous variables for the counterfactual
% scenarios.
create_counterfactual_scenarios



%% Model equilibrium
equilibrium_main



%% Volatility analysis
analysis
