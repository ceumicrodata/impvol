function data_management
%% Data management

global c

fe = c.fe;
fe_prices = c.fe_prices;
alpha_type = c.alpha_type;
i_base = c.i_base;
has_prices = c.has_prices;

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

io_values = dlmread([input_folder, 'oecd_io_values.csv'], ',', 1, 3);
total_output = dlmread([input_folder, 'oecd_total_output.csv'], ',', 1, 2);
output_shares = dlmread([input_folder, 'output_shares.csv'], ',', 1, 1);
intermediate_input_shares = dlmread([input_folder, 'intermediate_input_shares.csv'], ',', 1, 1);

trade_balance = dlmread([input_folder, 'trade_balance.csv'], ',', 1, 1);
trade_balance = 1000 * trade_balance; % convert to millions from billions

% trade_balance = dlmread([input_folder, 'trade_balance_new.csv'], ',', 1, 1);


%% Get key constants of data dimensions
n_countries = length(country_names);
n_sectors = size(beta_panel, 2);
n_years = size(beta_panel, 1) / n_countries;
i_services = n_sectors;

parameters.n_countries = n_countries;
parameters.n_sectors = n_sectors;
parameters.n_years = n_years;
parameters.i_services = i_services;


%% Reshape value added and betas
% Reshape value added to a more convenient 3 dimensional form with
% dimensions (country, sector, year).
va = reshape(reshape(va, n_years, n_countries * n_sectors)',...
    n_countries, n_sectors, n_years);


% To get sectoral betas take average across years and countries.
beta = mean(beta_panel)';
% beta_panel = reshape(reshape(beta_panel, n_years, n_countries * n_sectors)',...
%     n_countries, n_sectors, n_years);

%% Compute gammas
gammas = compute_gammas(io_values, total_output, output_shares, intermediate_input_shares);

% bar([1 - beta, sum(gammas, 1)'])

% normalize gammas
gammas = bsxfun(@times, gammas, (1 - beta') ./ sum(gammas, 1));

%% Compute alphas


alpha = compute_alphas(va, beta, gammas, weights);

gammas = repmat(gammas, [1, 1, n_years]);

% compute old alphas
% 
% if alpha_type == 0
%     % Constant alpha
%     va_sum_over_n_t = squeeze(sum(sum(va, 1), 3));
%     sectoral_va_share = va_sum_over_n_t / sum(va_sum_over_n_t);
%     alpha_old = (sectoral_va_share' ./ beta) / sum(sectoral_va_share' ./ beta);
%     alpha_old = repmat(alpha_old, [1 n_years]);
% else
%     % Time varying alpha
%     va_sum_over_n = squeeze(sum(va, 1));
%     sectoral_va_share = va_sum_over_n ./ ...
%                         repmat(sum(va_sum_over_n, 1), [n_sectors 1]);
%     sectoral_va_share_over_beta = sectoral_va_share ./ ...
%                                   repmat(beta, [1 n_years]);
%     
%     alpha_old = sectoral_va_share_over_beta ./ ...
%             repmat(sum(sectoral_va_share_over_beta, 1), [n_sectors 1]);
%       
%     % Smooth the time varying alphas using a band-pass filter
%     if alpha_type == 2
%         [alpha_trend, ~] = detrend_series(alpha_old, weights);
%         alpha_old = alpha_trend ./ repmat(sum(alpha_trend, 1), [n_sectors 1]);
%     end % if
% end % if

% alpha_const = repmat(mean(alpha, 2), [1, n_years]);
% 
% 
% alpha_linear = alpha;
% for j = 1:n_sectors
%     alpha_linear(j, :) = linspace(alpha(j, 1), alpha(j, n_years), n_years);
% end %j

% alpha = alpha_linear;
% alpha = alpha_const;

% alpha = repmat(mean(alpha_old, 2), [1, n_years]);
% gammas = mean(alpha, 2) * (1 - beta)';


% alpha = alpha_old;
% gammas = zeros(n_sectors, n_sectors, n_years);
% for t = 1:n_years
% gammas(:, :, t) = alpha(:, t) * (1 - beta)';
% end

%% Import sectoral Prices
% index of base country in sectoral price data
ii = sum(has_prices(1:i_base));
p_sectoral_data = csvread([input_folder, 'sectoral_price_index.csv']);
p_sectoral_base = p_sectoral_data((ii - 1) * n_years + 1 : ii * n_years, :);



%% Process imported matrices
% Normalize sectoral price index in base country 
p_sectoral_base = p_sectoral_base ./ repmat(p_sectoral_base(1, :), [n_years 1]);

% Compute aggregate price index in base country
p_base = prod((p_sectoral_base' ./ alpha) .^ alpha, 1)';

% Recompute PWT. Needed if US is not chosen as the base country
pwt = reshape(pwt, n_years, n_countries);
pwt = pwt ./ repmat(pwt(:, i_base), [1 n_countries]);

%% Convert import shares to expenditure shares based on source country.
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

B = bsxfun(@times, beta.^(-beta),  squeeze(prod(gammas.^(- gammas), 1)));

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
va_shares = va ./ permute(repmat(va_total, [1 1 n_sectors]), [2 3 1]);
[psi, ~] = detrend_series(va_shares, weights);


%% Calculate sectoral prices

p_sectoral = ...
    exp(squeeze(median(1/theta * log(bsxfun(@rdivide, d, d(i_base, :, :, :))) - ...
                log(bsxfun(@rdivide, kappa, kappa(i_base, :, :, :))), 2)) + ...
        permute(repmat(log(p_sectoral_base), [1, 1, n_countries]), [3, 2, 1]));
 
    
p_sectoral(:, i_services, :) = ...
    compute_p_services(pwt, p_sectoral, p_base, parameters);


% 
% %% check sectoral prices
% close all
% 
% for j = 1:24
%     figure()
% 
%     n = 2;
%     n1 = sum(has_prices(1:n));
%    
% 
%     a = p_sectoral_data(1 + (n1 - 1) * 36 : n1 * 36, j);
%     a = a / a(1);
% 
%     aa = squeeze(p_sectoral(n, j, :));
%     aa = aa / aa(1);
% 
%     plot([a, aa])
% 
% end

%% get Z

z = zeros(n_countries, n_sectors, n_years);


for n = 1:n_countries
    for t = 1:n_years
        z(n, :, t) = theta * (beta' .* (log(va(n, :, t)) - log(psi(n, :, t))) + log(p_sectoral(n, :, t)) * gammas(:, :, t) + log(xi * B(:, t)')) + ...
                     mean(squeeze(log(d(:, n, :, t)) - theta * log(kappa(:, n, :, t))) - theta * log(p_sectoral(:, :, t)), 1);
    end
end


z = exp(z);


%% Calculate productivity shocks

% 
% zeta = compute_zeta(d, va, psi, pwt, p_sectoral, parameters);
% z_traded = exp(squeeze(mean(zeta - theta * permute(repmat(log(p_sectoral(:, 1:(n_sectors - 1), :)), [1 1 1 n_countries]), [1 4 2 3]), 2)));
% 
% z = zeros(n_countries, n_sectors, n_years);
% z(:, 1:(n_sectors - 1), :) = z_traded;

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
% baseline.K = K;
baseline.alpha = alpha;
baseline.beta = beta;
baseline.kappa = kappa;
baseline.theta = theta;
baseline.xi = xi;
baseline.gammas = gammas;


deflator = pwt .* repmat(p_base, [1, n_countries]);
% baseline.trade_balance = bsxfun(@rdivide, trade_balance, p_base');
% baseline.trade_balance = trade_balance ./ deflator';
baseline.trade_balance = trade_balance;

baseline.L = L;
baseline.z = z;
baseline.psi = psi;

save([output_folder, 'data_theta_', num2str(theta), '.mat'], 'baseline');



%% Compute and save Real GDP and volatility in data

real_gdp_total = (va_total ./ deflator)';
real_gdp_sectoral = va ./ permute(repmat(deflator', [1 1 n_sectors]), [1 3 2]);

[~, data_cycle_total] = detrend_series(log(real_gdp_total), weights);
data_volatility_total = var(data_cycle_total, 0, 2);

output_folder = c.results_folder;
save([output_folder, 'data_rgdp_and_volatility.mat'], 'country_names', 'real_gdp_sectoral',...
     'real_gdp_total', 'data_volatility_total', 'deflator', 'pwt', 'p_base',...
     'va_total', 'va', 'p_sectoral_base', 'p_sectoral_data', 'd');
end