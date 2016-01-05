%% Main file to replicate the impvol project
%% Preliminaries
clear all
clc
format compact
addpath(genpath(pwd))

%tic

global c

%% Initialize program parameters
c = initialize4;



%% Data management
% Import and process data:
% Compute parameters and exogenous variables of the model.
data_management



%% Create counterfactual scenarios
% Compute alternative parameters and exogenous variables for the counterfactual
% scenarios.
create_counterfactual_scenarios_1972_actual



%% Model equilibrium
equilibrium_main


%% Analysis
analysis

%toc