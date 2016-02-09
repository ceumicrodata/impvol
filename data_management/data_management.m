function data_management
%% Data management

global c

fe = c.fe;
fe_prices = c.fe_prices;
alpha_type = c.alpha_type;
i_base = c.i_base;
have_prices = c.have_prices;

theta = c.theta;
eta = c.eta;
numerical_zero = c.numerical_zero;
weights = c.filter_weights;

input_folder = c.data_folder_original;
output_folder = c.data_folder_algorithm_input;

parameters.theta = theta;
parameters.numerical_zero = numerical_zero;
parameters.i_base = i_base;



%% Import files from data folder.
country_names = importdata([input_folder, 'country_name.txt']);
beta_panel = dlmread([input_folder, 'beta_panel.txt'], '\t', 1, 2);

pwt = dlmread([input_folder,...
    'aggregate_price_relative_to_US.csv'], ',', 0, 2);

va = dlmread([input_folder, 'sectoral_value_added.txt'], '\t', 1, 2);
import_shares = dlmread([input_folder, 'import_share.txt'], '\t');

% gammas = dlmread([input_folder, 'io_linkages_parameters.txt']);


%% Get key constants of data dimensions
n_countries = length(country_names);
n_sectors = size(beta_panel, 2);
n_years = size(beta_panel, 1) / n_countries;
i_services = n_sectors;

parameters.n_countries = n_countries;
parameters.n_sectors = n_sectors;
parameters.n_years = n_years;
parameters.i_services = i_services;



%% Compute alphas
% Reshape value added to a more convenient 3 dimensional form with
% dimensions (country, sector, year).
va = reshape(reshape(va, n_years, n_countries * n_sectors)',...
    n_countries, n_sectors, n_years);

% To get sectoral betas take average across years and countries.
beta = mean(beta_panel)';

if alpha_type == 0
    % Constant alpha
    va_sum_over_n_t = squeeze(sum(sum(va, 1), 3));
    sectoral_va_share = va_sum_over_n_t / sum(va_sum_over_n_t);
    alpha = (sectoral_va_share' ./ beta) / sum(sectoral_va_share' ./ beta);
    alpha = repmat(alpha, [1 n_years]);
else
    % Time varying alpha
    va_sum_over_n = squeeze(sum(va, 1));
    sectoral_va_share = va_sum_over_n ./ ...
                        repmat(sum(va_sum_over_n, 1), [n_sectors 1]);
    sectoral_va_share_over_beta = sectoral_va_share ./ ...
                                  repmat(beta, [1 n_years]);
    
    alpha = sectoral_va_share_over_beta ./ ...
            repmat(sum(sectoral_va_share_over_beta, 1), [n_sectors 1]);
      
    % Smooth the time varying alphas using a band-pass filter
    if alpha_type == 2
        [alpha_trend, ~] = detrend_series(alpha, weights);
        alpha = alpha_trend ./ repmat(sum(alpha_trend, 1), [n_sectors 1]);
    end % if
end % if



%% Import sectoral Prices
% index of base country in sectoral price data
ii = sum(have_prices(1:i_base));
p_sectoral_data = csvread([input_folder, 'sectoral_price_index.csv']);
p_sectoral_base = p_sectoral_data((ii - 1) * n_years + 1 : ii * n_years, :);



%% Process imported matrices
% Normalize sectoral price index in base country 
p_sectoral_base = p_sectoral_base ./ repmat(p_sectoral_base(1, :), [n_years 1]);

% Compute aggregate price index in base country
p_base = prod((p_sectoral_base' ./ alpha) .^ alpha)';

% Recompute PWT. Needed if US is not chosen as the base country
pwt = reshape(pwt, n_years, n_countries);
pwt = pwt ./ repmat(pwt(:, i_base), [1 n_countries]);

% Convert import shares to expenditure shares based on source country.
d = zeros(n_countries, n_countries, n_sectors, n_years);
n_lines_in_import_shares = length(import_shares);

first_year = import_shares(1, 1);
for i_line = 1:n_lines_in_import_shares
    year = import_shares(i_line, 1) + 1 - first_year;
    importer = import_shares(i_line, 2);
    exporter = import_shares(i_line, 3);
    
    % Clip negative import shares from below at 0.
    d(importer, exporter, 1:(n_sectors - 1), year) = ...
        max(import_shares(i_line, 4:26), 0);
    
    % Set services shares to 0 as they are not traded.
    d(importer, exporter, i_services, year) = 0;
end

% There are a number of manipulations to the shares
d = correct_expenditure_shares(d, parameters);


%% Calculate new matrices of constants.
xi = gamma((theta + 1 - eta) / theta);


% experimental gamma
rng(0);
gammas = rand(n_sectors, n_sectors);
gammas = bsxfun(@times, gammas, (1 - beta') ./ sum(gammas));


B = beta.^(-beta) .* prod(gammas.^(- gammas), 1)';

K = zeros(n_years, 1);
for t = 1:n_years
    K(t) = prod((xi * B ./ alpha(:, t)).^(alpha(:, t)))^theta;
end % t

va_total = sum(va, 2);
va_total = squeeze(va_total)';

kappa = compute_trade_cost(d, parameters);


% collect parameters
parameters.alpha = alpha;
parameters.B = B;
parameters.beta = beta;
parameters.xi = xi;
parameters.kappa = kappa;
parameters.gammas = gammas;


%% Get expectations of value added shares (psi)
va_shares = va ./ permute(repmat(va_total', [1 1 n_sectors]), [1 3 2]);
[psi, ~] = detrend_series(va_shares, weights);



%% Calculate the log difference in variables
p_sectoral_base_change = diff(log(p_sectoral_base));
pwt_change = diff(log(pwt));
p_base_change = diff(log(p_base));

% va_totalChange = diff(log(va_total));
va_change = diff(log(va), 1, 3);
psi_change = diff(log(psi), 1, 3);

kappa_to_base_change = diff(log(squeeze(kappa(i_base, :, :, :))), 1, 3);
d_to_base_change = diff(log(squeeze(d(i_base, :, :, :))), 1, 3);



%% Calculate sectoral prices

p_sectoral = ...
    exp(squeeze(mean(1/theta * bsxfun(@minus, log(d), log(d(i_base, :, :, :))) - ...
                bsxfun(@minus, log(kappa), log(kappa(i_base, :, :, :))), 2)) + ...
        permute(repmat(log(p_sectoral_base), [1, 1, n_countries]), [3, 2, 1]));

p_sectoral(:, i_services, :) = ...
    compute_p_services(pwt, p_sectoral, p_base, parameters);


%% Calculate productivity shocks



%%

zeta = compute_zeta(d, va, psi, pwt, p_sectoral, parameters);
z_traded = exp(squeeze(mean(zeta - theta * permute(repmat(log(p_sectoral(:, 1:(n_sectors - 1), :)), [1 1 1 n_countries]), [1 4 2 3]), 2)));

z = zeros(n_countries, n_sectors, n_years);
z(:, 1:(n_sectors - 1), :) = z_traded;

% Compute z of services in 1972 in all countries
z(:, i_services, :) = ...
    compute_z_services(va, psi, p_sectoral, parameters);


%% Compute Equipped Labor - L
% This is alpha_j * (37) summed over j and then applying (38).
full_power = permute(repmat(alpha ./ repmat(beta * theta, [1 n_years]),...
                            [1, 1, n_countries]),...
                     [3 1 2]);
L = prod(z.^full_power, 2);
L = squeeze(L);



%% Save computed variables and parameters
baseline = struct;
baseline.scenario = 'Baseline';

baseline.country_names = country_names;
baseline.B = B;
baseline.K = K;
baseline.alpha = alpha;
baseline.beta = beta;
baseline.kappa = kappa;
baseline.theta = theta;
baseline.xi = xi;
baseline.gammas = gammas;

baseline.L = L;
baseline.z = z;
baseline.psi = psi;

save([output_folder, 'data_theta_', num2str(theta), '.mat'], 'baseline');



%% Compute and save Real GDP and volatility in data
% deflator = pwt .* repmat(pUs * pUs1972, [1, n_countries]);
deflator = pwt .* repmat(p_base, [1, n_countries]);

real_gdp_total = (va_total ./ deflator)';
real_dgp_sectoral = va ./ permute(repmat(deflator', [1 1 n_sectors]), [1 3 2]);

[~, data_cycle_total] = detrend_series(log(real_gdp_total), weights);
data_volatility_total = var(data_cycle_total, 0, 2);

output_folder = c.results_folder;
save([output_folder, 'data_rgdp_and_volatility.mat'], 'real_dgp_sectoral',...
     'real_gdp_total', 'data_volatility_total', 'deflator', 'pwt', 'p_base',...
     'va_total', 'va', 'p_sectoral_base', 'p_sectoral_data', 'd');
end