function c = initialize_ubt

%% Specify run type
% 0: no counterfactuals, 1: do counterfactuals
c.include_counterfactuals = 1; 

c.io_links = 0;

c.ubt = 1;

% 0: const, 1: general(time-varying), 2: smooth
c.alpha_type = 2; 

% 0: Use structural model equations to compute the productivities (z)
% 1: Use the fixed effects approach to compute the productivities (z)
c.fe = 1; 

% 0: Use structural model equations to compute sectoral prices
% 1: Use sectoral prices from fixed effects decomposition
c.fe_prices = 1; 
assert(c.fe >= c.fe_prices,...
       'Can''t have fe prices without fe decomposition')



%% Set main parameters
% Set main model parameters
c.theta = 4;
c.eta = 2;

% Specify base country. (Has to have sectoral prices available)
c.i_base = 25;

% Specify countries that have sectoral price index
c.has_prices = ...
    logical([1 1 1 0 0 0 1 1 1 1 1 0 1 1 1 0 1 0 1 0 1 1 1 1 1])';
assert(c.has_prices(c.i_base) == 1, 'No sectoral prices for base country.')

% Set folder locations
c.data_folder_original = 'data/raw_imputed/';
c.data_folder_algorithm_input = ['data/algorithm_input_io_', num2str(c.io_links), '_ubt_', num2str(c.ubt), '/'];
c.results_folder = ['results_io_', num2str(c.io_links), '_ubt_', num2str(c.ubt), '/'];

mkdir(c.data_folder_algorithm_input)
mkdir(c.results_folder)

% Set amount of info printed on screen
% 0: no output
% 1: labor loop iterations (show diff)
% 2: 1 + wage loop convergence info for all periods (show number of steps)
% 3: 2 + wage loop iteration info (diff for every nth iteration) 
% 4: 3 + price loop convergence info (number of steps + mean aggr. price)
c.verbosity = 3;
c.outer_print_every = 1;

% specify how frequency of iteration display in wage loop
c.middle_print_every = 25;

% Band pass filter weights
c.filter_weights = [0.774074394803123;...
                  -0.201004684236153;...
                  -0.135080548288772;...
                  -0.050951964876636];
              
% Or put [0] for log differencing
% c.filter_weights = [0];


c.numerical_zero = 1e-12;

% Dampening parameters for the equilibrium search loops
c.dampening_labor_loop = 0.75;
c.dampening_wage_loop = 0.4:-0.05:0.2;
c.dampening_price_loop = 0.7:-0.05:0.1;

% Technical values for the loops
c.dif = 1e9; % a big number

% Technical values for the outer loop
c.outer_tol = 1e-2; % set the convergence tolerance
c.outer_maxiter = 50; % limit the maximum number of iterations

% Technical values for the middle loop
c.middle_tol = 1e-3; % set the convergence tolerance
c.middle_maxiter = 4e2; % limit the maximum number of iterations

% Technical values for the inner loop
c.inner_tol = 1e-3; % set the convergence tolerance
c.inner_maxiter = 2e2; % limit the maximum number of iterations
end