function c = init_globals(model)

%% Load values from specification
load(['model_specifications/', model, '.mat'])
c = spec;

% Specify base country. (Has to have sectoral prices available)
c.i_base = 25;

% Specify countries that have sectoral price index
c.has_prices = ...
    logical([1 1 1 0 0 0 1 1 1 1 1 0 1 1 1 0 1 0 1 0 1 1 1 1 1])';
assert(c.has_prices(c.i_base) == 1, 'No sectoral prices for base country.')


% Set folder locations
% example target folder model00

% fname = [model_id,...
%         '_iol_', num2str(c.iol),...
%         '_bt_', num2str(c.bt),...
%         '_lac_', num2str(c.lac),...
%         '_th_', num2str(c.th),...
%         '_et_', num2str(c.et),...
%         '_tc_', num2str(c.tc),...
%         '_sh_', num2str(c.sh)];

c.data_folder_original = 'data/raw_imputed/';
c.model_folder = ['models/', model, '/'];

mkdir(c.model_folder)


% Set amount of info printed on screen
% 0: no output
% 1: labor loop iterations (show diff)
% 2: 1 + wage loop convergence info for all periods (show number of steps)
% 3: 2 + wage loop iteration info (diff for every nth iteration) 
% 4: 3 + price loop convergence info (number of steps + mean aggr. price)
c.verbosity = 2;
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
c.dampening_labor_loop = 1;
c.dampening_wage_loop = 0.1;
c.dampening_price_loop = 0.25;

% Technical values for the loops
c.dif = 1e9; % a big number - really unnecessary...

% Technical values for the outer loop
c.outer_tol = 1e-3; % set the convergence tolerance
c.outer_maxiter = 1e2; % limit the maximum number of iterations

% Technical values for the middle loop
c.middle_tol = 1e-4; % set the convergence tolerance
c.middle_maxiter = 2e3; % limit the maximum number of iterations

% Technical values for the inner loop
c.inner_tol = 1e-4; % set the convergence tolerance
c.inner_maxiter = 2e3; % limit the maximum number of iterations


% save model specification in model folder
save([c.model_folder, 'specification.mat'], '-struct', 'c');

end
