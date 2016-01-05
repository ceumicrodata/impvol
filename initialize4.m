function c = initialize4

%% Specify run type
% 0: no counterfactuals, 1: do counterfactuals
c.includeCounterfactuals = 1; 

% 0: const, 1: general(time-varying), 2: smooth
c.alphaType = 2; 

% 0: Use structural model equations to compute the productivities (Z)
% 1: Use the fixed effects approach to compute the productivities (Z)
c.FE = 1; 

% 0: Use structural model equations to compute sectoral prices
% 1: Use sectoral prices from FE decomposition
c.FE_prices = 1; 
assert(c.FE >= c.FE_prices,...
       'Can''t have FE prices without FE decomposition')



%% Set main parameters
% Set main model parameters
c.theta = 4;
c.eta = 2;

% Specify base country. (Has to have sectoral prices available)
c.iBase = 25;

% Specify countries that have sectoral price index
c.havePrices = ...
    logical([1 1 1 0 0 0 1 1 1 1 1 0 1 1 1 0 1 0 1 0 1 1 1 1 1])';
assert(c.havePrices(c.iBase) == 1, 'No sectoral prices for base country.')

% Set folder locations
c.dataFolderOriginal = 'data/raw_imputed/';
c.dataFolderAlgorithmInput = 'data/algorithm_input/';
c.resultsFolder = 'results/';

% Set amount of info printed on screen
% 0: no output
% 1: labor loop iterations (show diff)
% 2: 1 + wage loop convergence info for all periods (show number of steps)
% 3: 2 + wage loop iteration info (diff for every nth iteration) 
% 4: 3 + price loop convergence info (number of steps + mean aggr. price)
c.verbosity = 1;
c.outerPrintEvery = 1;

% specify how frequency of iteration display in wage loop
c.middlePrintEvery = 20;

% Band pass filter weights
c.filterWeights = [0.774074394803123;...
                  -0.201004684236153;...
                  -0.135080548288772;...
                  -0.050951964876636];
              
% Or put [0] for log differencing
% c.filterWeights = [0];


c.numericalZero = 1e-10;

% Dampening parameters for the equilibrium search loops
c.dampeningLaborLoop = 0.75;
c.dampeningWageLoop = 0.55:-0.05:0.2;
c.dampeningPriceLoop = 0.7:-0.05:0.1;

% Technical values for the outer loop
c.outer_dif = 1e9; % a big number
c.outer_tol = 1e-3; % set the convergence tolerance
c.outer_maxiter = 250; % limit the maximum number of iterations

% Technical values for the middle loop
c.middle_dif = 1e9; % a big number
c.middle_tol = 1e-3; % set the convergence tolerance
c.middle_maxiter = 1e3; % limit the maximum number of iterations

% Technical values for the inner loop
c.inner_dif = 1e9; % a big number
c.inner_tol = 1e-3; % set the convergence tolerance
c.inner_maxiter = 5e3; % limit the maximum number of iterations
end