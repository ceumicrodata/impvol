clear all
%% Initialize program parameters
for io_links = 0:1
	for unbalanced_trade = 1-io_links:1


		%% Main file to replicate the impvol project
		%% Preliminaries
		%%clear all
		clc
		format compact
		addpath(genpath(pwd))

		theta = 4;

		%tic

		global c

		c = initialize_parameters(theta, io_links, unbalanced_trade);



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
	end
end
