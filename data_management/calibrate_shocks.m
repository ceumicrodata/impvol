function calibrate_shocks
%% Data management

global c

i_base = c.i_base;
io_links = c.iol;
bt = c.bt;

theta = c.th;
eta = c.et;
rho = c.rh;
numerical_zero = c.numerical_zero;
weights = c.filter_weights;

input_folder = c.data_folder_original;

parameters.rho = rho;
parameters.theta = theta;
parameters.numerical_zero = numerical_zero;
parameters.i_base = i_base;
parameters.has_prices = c.has_prices;

%% Import files from data folder.
country_names = importdata([input_folder, 'country_name.txt']);
%sector_names = importdata([input_folder, 'beta_panel.txt']);
%sector_names = sector_names.textdata(1,3:end);
beta_panel = wrapper(dlmread([input_folder, 'beta_panel.txt'], '\t', 1, 2));

% va = wrapper(dlmread([input_folder, 'sectoral_value_added.txt'], '\t', 1, 2));
va = wrapper(dlmread([input_folder, 'sectoral_value_added.csv'], ',', 1, 2));

import_shares = wrapper(dlmread([input_folder, 'import_share.txt'], '\t'));

io_values = wrapper(dlmread([input_folder, 'oecd_io_values.csv'], ',', 1, 3));
total_output = wrapper(dlmread([input_folder, 'oecd_total_output.csv'], ',', 1, 2));
output_shares = wrapper(dlmread([input_folder, 'output_shares.csv'], ',', 1, 1));
intermediate_input_shares = wrapper(dlmread([input_folder, 'intermediate_input_shares.csv'], ',', 1, 1));

trade_balance = wrapper(dlmread([input_folder, 'trade_balance_new.csv'], ',', 1, 1));
trade_balance = bsxfun(@minus, trade_balance, mean(trade_balance));

if bt == 1
    trade_balance = zeros(size(trade_balance));
end %if

p_sectoral_data = csvread([input_folder, 'sectoral_price_index.csv']);
pwt = wrapper(dlmread([input_folder,...
    'aggregate_price_relative_to_US.csv'], ',', 0, 2));

%% Get data dimensions
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
    n_countries, n_sectors, n_years); % (country, sector, year)

va_share_orig = va./repmat(sum(va,2),[1,n_sectors,1]);

va_long = long(va,parameters); % (country×sector,year)

% Get expectations of value added shares (psi)
va_total = squeeze(sum(va, 2))';

va_shares = va ./ permute(repmat(va_total, [1 1 n_sectors]), [2 3 1]);
[psi, ~] = detrend_series(va_shares, weights);

% To get sectoral betas take average across years and countries.
beta_replace = mean(beta_panel)'; %(sector,1) - country and year average
% Here we calculate the betas as an average over countries - this could be
% done prettier
beta = zeros(n_countries,n_sectors,n_years);
for t = 1:n_years;
    beta(:,:,t) = beta_panel(t:n_years:(n_years - 1)*n_countries + t,:);
end
beta = squeeze(mean(beta,1)); %(sector,year) - country average


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
d = correct_expenditure_shares(d, parameters); % (imp_country, exp_country, sector, year)

D = zeros(n_countries*n_sectors,n_countries*n_sectors,n_years);
for sector = 1:n_sectors
    exporter = sector:n_sectors:(n_countries-1)*n_sectors + sector;
    for imp = 1:n_countries
        for t = 1:n_years
            D(n_sectors*(imp - 1) + sector, exporter,t) = d(imp,:,sector,t); % (imp_country×sector, exp_country×sector, year)
        end
    end
    clear exp
end

%% IO links section
if io_links == 1
    gamma_x = compute_gammas(io_values, total_output, output_shares, intermediate_input_shares);
    gammas = zeros(n_sectors, n_sectors, n_years);
    for t = 1:n_years
        gammas(:,:,t) = bsxfun(@times, gamma_x, (1 - beta(:,t)') ./ sum(gamma_x, 1));
    end
    
    final_expenditure_share = compute_final_expenditure_share(D, va_long, gammas, beta, parameters, c);
else % this branch is broken
    va_sum_over_n = squeeze(sum(va, 1));
    sectoral_va_share = va_sum_over_n ./ ...
        repmat(sum(va_sum_over_n, 1), [n_sectors 1]);
    sectoral_va_share_over_beta = sectoral_va_share ./ ...
        repmat(beta_replace, [1 n_years]);
    alpha_old = sectoral_va_share_over_beta ./ ...
        repmat(sum(sectoral_va_share_over_beta, 1), [n_sectors 1]);
    
    [alpha_trend, ~] = detrend_series(alpha_old, weights);
    alpha_old = alpha_trend ./ repmat(sum(alpha_trend, 1), [n_sectors 1]);
    
    alpha_replace = alpha_old;
    gammas_replace = zeros(n_sectors, n_sectors, n_years);
    for t = 1:n_years
        gammas_replace(:, :, t) = alpha_replace(:, t) * (1 - beta_replace)';
    end
end

%% Calculate new matrices of constants.
theta = c.th;
eta = c.et;
xi = gamma((theta + 1 - eta) / theta);

B = beta.^(-beta) .* squeeze(prod(gammas.^(-gammas), 1));
assert_all(size(B)==[n_sectors, n_years]);

kappa = compute_trade_cost(d, parameters);
assert_all(kappa>=0);
assert_all(kappa<=1);

% collect parameters
parameters.final_expenditure_share = final_expenditure_share;
parameters.kappa = kappa;

%% Import sectoral Prices

% Recompute PWT. Needed if US is not chosen as the base country
pwt = reshape(pwt, n_years, n_countries);
pwt = pwt ./ repmat(pwt(:, i_base), [1 n_countries]);

% Calculate sectoral prices
[p_sectoral, p_base, p_sectoral_base] = compute_sectoral_p(p_sectoral_data, pwt, d, parameters);

%% Recover nu from alpha using rho, nu = alpha * p_sectoral^(rho - 1) from alpha = normalization_term * nu * p_sectoral^(1 - rho). alpha must sum to unity over sectors
nu = compute_nus(p_sectoral, parameters);


% % collect parameters
% parameters.final_expenditure_share = final_expenditure_share;
% parameters.nu = nu;
% parameters.B = B;
% parameters.beta = beta;
% parameters.xi = xi;
% parameters.kappa = kappa;
% parameters.gammas = gammas;

%% get Z

z = zeros(n_countries, n_sectors, n_years);

for n = 1:n_countries
    for t = 1:n_years
        z(n, :, t) = theta * log(xi * B(:, t)') + ...
                     theta * beta(:,t)' .* (log(va(n, :, t)) - log(psi(n, :, t))) + ...
                     mean(squeeze(log(d(:, n, :, t)) - theta * log(kappa(:, n, :, t))), 1) - ...
                     theta * mean(log(p_sectoral(:, :, t)), 1) + ...
                     theta * log(p_sectoral(n, :, t)) * gammas(:, :, t);
    end % for t
end % for n


z = exp(z);


%% Calculate productivity shocks

z_services = zeros(n_countries, n_years);

for n = 1:n_countries
    for t = 1:n_years
        z_services(n, t) = ...
            xi^theta * ...
            B(i_services, t)^theta * ...
            (va(n, i_services, t) / psi(n, i_services, t))^...
                                            (theta * beta(i_services,t)) * ...
            prod(p_sectoral(n, :, t).^(gammas(:, i_services, t)'))^theta * ...
            p_sectoral(n, i_services, t)^(-theta);
    end % t
end % n

z(:, i_services, :) = z_services;

%% Compute Equipped Labor - L
% Verify with scaling, it should not matter.
L = ones(n_countries, n_years);

%% Save computed variables and parameters


baseline = struct;
load([c.model_folder, 'specification.mat'], 'model');
baseline.scenario = model;

baseline.country_names = country_names;
baseline.B = B;
% baseline.K = K;
baseline.final_expenditure_share = final_expenditure_share;
baseline.nu = nu;
baseline.beta = beta;
baseline.kappa = kappa;
baseline.theta = theta;
baseline.xi = xi;
baseline.rho = rho;
baseline.gammas = gammas;

baseline.trade_balance = trade_balance;

baseline.L = L;
baseline.z = z;
baseline.psi = psi;

save([c.model_folder, 'alg_inputs.mat'], 'baseline');



%% Compute and save Real GDP and volatility in data

deflator = pwt .* repmat(p_base, [1, n_countries]);
% baseline.trade_balance = bsxfun(@rdivide, trade_balance, p_base');
% baseline.trade_balance = trade_balance ./ deflator';

real_gdp_total = (va_total ./ deflator)';
real_gdp_sectoral = va ./ permute(repmat(deflator', [1 1 n_sectors]), [1 3 2]);

[~, data_cycle_total] = detrend_series(log(real_gdp_total), weights);
data_volatility_total = var(data_cycle_total, 0, 2);

save([c.model_folder, 'data_rgdp_and_volatility.mat'], 'country_names',...
     'data_volatility_total', 'deflator', 'pwt', 'p_base',...
     'va_total', 'va', 'p_sectoral_base', 'p_sectoral_data', 'p_sectoral', 'd');
end