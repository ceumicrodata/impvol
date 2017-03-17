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

parameters.theta = theta;
parameters.numerical_zero = numerical_zero;
parameters.i_base = i_base;



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
%Inverse reshape: va2 = permute(reshape(VA,n_sectors, n_countries, n_years),[2,1,3]);
%for imp = 1:n_countries
%    for t = 1:n_years
%        VA(n_sectors*(imp - 1) + 1:n_sectors*imp, n_sectors*(imp - 1) + 1:n_sectors*imp, t) = repmat(va(imp,:,t),[n_sectors,1]); % Am I sure about the order of the dimensions?
%    end
%end

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

% Identity matrix
% D = eye(n_countries*n_sectors);
% D = repmat(D,[1,1,n_years]);

%% IO links section
if io_links == 1
    gamma = compute_gammas(io_values, total_output, output_shares, intermediate_input_shares);
    gammas = zeros(n_sectors, n_sectors, n_years);
    for t = 1:n_years
        gammas(:,:,t) = bsxfun(@times, gamma, (1 - beta(:,t)') ./ sum(gamma, 1));
    end
    clear gamma;
    
    final_expenditure_share = compute_final_expenditure_share(D, va_long, gammas, beta, parameters, c);
    
    %final_expenditure_share_long = long(final_expenditure_share, parameters);    
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

%% Get expectations of value added shares (psi)
va_total = squeeze(sum(va, 2))';

va_shares = va ./ permute(repmat(va_total, [1 1 n_sectors]), [2 3 1]);
[psi, ~] = detrend_series(va_shares, weights);

%% Calculate new matrices of constants.
theta = c.th;
eta = c.et;
xi = gamma((theta + 1 - eta) / theta);

%This is a constant in case of Cobb-Douglas production function - original 
B_replace = bsxfun(@times, beta_replace.^(-beta_replace),  squeeze(prod(gammas_replace.^(-gammas_replace), 1)));
%New B, shouldn't it be something else?
B = bsxfun(@times, beta.^(-beta),  squeeze(prod(gammas.^(-gammas), 1)));

kappa = compute_trade_cost(d, parameters);

% collect parameters
parameters.alpha_r = alpha_replace;
parameters.alpha_t = final_exp_share_tr_matr;
parameters.B_r = B_replace;
parameters.B_t = B;
parameters.beta_r = beta_replace;
parameters.beta_t = beta;
parameters.xi = xi;
parameters.kappa = kappa;
parameters.gammas_r = gammas_replace;
parameters.gammas_t = gammas;

%% Import sectoral Prices

% Recompute PWT. Needed if US is not chosen as the base country
pwt = reshape(pwt, n_years, n_countries);
pwt = pwt ./ repmat(pwt(:, i_base), [1 n_countries]);

% Calculate sectoral prices
[p_sectoral, p_base, p_sectoral_base] = compute_sectoral_p(p_sectoral_data, pwt, d, c, parameters);
% % index of base country in sectoral price data
% ii = sum(has_prices(1:i_base));
% p_sectoral_data = csvread([input_folder, 'sectoral_price_index.csv']);
% p_sectoral_base = p_sectoral_data((ii - 1) * n_years + 1 : ii * n_years, :);
% 
% % aa = cumsum(has_prices);
% % 
% % p_sectoral_data2 = zeros(n_countries, n_sectors, n_years);
% % for n = 1:n_countries
% %     if has_prices(n) == 1
% %         p_sectoral_data2(n, :, :) =  p_sectoral_data((aa(n) - 1) * n_years + 1 : aa(n) * n_years, : )';
% % %         p_sectoral_data2(n, :, :) = bsxfun(@rdivide, p_sectoral_data2(n, :, :), p_sectoral_data2(n, :, 1));  
% %     end % if
% % end % for n
% 
% 
% % p_sectoral_data2 = bsxfun(@rdivide, p_sectoral_data2, p_sectoral_data2(25, :, 1));

%% Process imported matrices
% Normalize sectoral price index in base country 
% p_sectoral_base = p_sectoral_base ./ repmat(p_sectoral_base(1, :), [n_years 1]);
% 
% % Compute aggregate price index in base country
% p_base = prod((p_sectoral_base' ./ alpha_replace) .^ alpha_replace, 1)';
% 
% % Recompute PWT. Needed if US is not chosen as the base country
% pwt = reshape(pwt, n_years, n_countries);
% pwt = pwt ./ repmat(pwt(:, i_base), [1 n_countries]);
% 
% %% Calculate new matrices of constants.
% xi = gamma((theta + 1 - eta) / theta);
% 
% B = bsxfun(@times, beta.^(-beta),  squeeze(prod(gammas.^(-gammas), 1)));
% 
% va_total = sum(va, 2);
% va_total = squeeze(va_total)';
% 
% kappa = compute_trade_cost(d, parameters);
% 
% % collect parameters
% parameters.alpha = alpha_replace;
% parameters.B = B;
% parameters.beta = beta;
% parameters.xi = xi;
% parameters.kappa = kappa;
% parameters.gammas = gammas;
%
% %% Calculate sectoral prices
% 
% p_sectoral = ...
%     exp(squeeze(mean(1/theta * log(bsxfun(@rdivide, d, d(i_base, :, :, :))) - ...
%                 log(bsxfun(@rdivide, kappa, kappa(i_base, :, :, :))), 2)) + ...
%         permute(repmat(log(p_sectoral_base), [1, 1, n_countries]), [3, 2, 1]));
%     
% p_sectoral(:, i_services, :) = ...
%     compute_p_services(pwt, p_sectoral, p_base, parameters);
% 
% 
% p_sectoral(:,i_services, :) = 1;
% % nu(services) = 1 for all n,t,
% % other sectoral prices normalized to 1 initially
% 
% % 
% % %% check sectoral prices
% % close all
% % 
% % for j = 1:24
% %     figure()
% % 
% %     n = 2;
% %     n1 = sum(has_prices(1:n));
% %    
% % 
% %     a = p_sectoral_data(1 + (n1 - 1) * 36 : n1 * 36, j);
% %     a = a / a(1);
% % 
% %     aa = squeeze(p_sectoral(n, j, :));
% %     aa = aa / aa(1);
% % 
% %     plot([a, aa])
% % 
% % end


%% Recover nu from alpha using rho, nu = alpha * p_sectoral^(rho - 1) from alpha = normalization_term * nu * p_sectoral^(1 - rho). alpha must sum to unity over sectors
nu =copmute_nus(final_exp_share, p_sectoral, parameters, c);

%%
% 
% n = 22;
% j = 21;
% 
% plot([squeeze(p_sectoral_data2(n, j, :)), squeeze(p_sectoral(n, j, :))])



%% get Z

z = zeros(n_countries, n_sectors, n_years);
% factor1 = zeros(size(z));
% factor2 = zeros(size(z));
% factor3 = zeros(size(z));
% factor4 = zeros(size(z));
% factor5 = zeros(size(z));

%xi = parameters.xi;
%B = parameters.B;
%kappa = parameters.kappa;

for n = 1:n_countries
    for t = 1:n_years
        z(n, :, t) = theta * log(xi * B_replace(:, t)') + ...
                     theta * beta_replace' .* (log(va(n, :, t)) - log(psi(n, :, t))) + ...
                     mean(squeeze(log(d(:, n, :, t)) - theta * log(kappa(:, n, :, t))), 1) - ...
                     theta * mean(log(p_sectoral(:, :, t)), 1) + ...
                     theta * log(p_sectoral(n, :, t)) * gammas_replace(:, :, t);
                 
%         for j = 1:n_sectors
%             factor1(n, j, t) = theta * log(xi * B(j, t));
%             factor2(n, j, t) = theta * beta(j) * (log(va(n, j, t)) - log(psi(n, j, t)));
%             factor3(n, j, t) = mean(log(d(:, n, j, t)) - theta * log(kappa(:, n, j, t)));
%             factor4(n, j, t) = - theta * mean(log(p_sectoral(:, j, t)));
%             factor5(n, j, t) = theta * log(p_sectoral(n, :, t)) * gammas(:, j, t);
%         end % for j
    end % for t
end % for n

for n = 1:n_countries
    for t = 1:n_years
        z(n, :, t) = theta * log(xi * B(:, t)') + ...
                     theta * beta(:,t)' .* (log(va(n, :, t)) - log(psi(n, :, t))) + ...
                     mean(squeeze(log(d(:, n, :, t)) - theta * log(kappa(:, n, :, t))), 1) - ...
                     theta * mean(log(p_sectoral(:, :, t)), 1) + ...
                     theta * log(p_sectoral(n, :, t)) * Gammas(:, :, t);
                 
%         for j = 1:n_sectors
%             factor1(n, j, t) = theta * log(xi * B(j, t));
%             factor2(n, j, t) = theta * beta(j) * (log(va(n, j, t)) - log(psi(n, j, t)));
%             factor3(n, j, t) = mean(log(d(:, n, j, t)) - theta * log(kappa(:, n, j, t)));
%             factor4(n, j, t) = - theta * mean(log(p_sectoral(:, j, t)));
%             factor5(n, j, t) = theta * log(p_sectoral(n, :, t)) * gammas(:, j, t);
%         end % for j
    end % for t
end % for n


z = exp(z);


%% Calculate productivity shocks

% 
% zeta = compute_zeta(d, va, psi, pwt, p_sectoral, parameters);
% z_traded = exp(squeeze(mean(zeta - theta * permute(repmat(log(p_sectoral(:, 1:(n_sectors - 1), :)), [1 1 1 n_countries]), [1 4 2 3]), 2)));
% 
% z = zeros(n_countries, n_sectors, n_years);
% z(:, 1:(n_sectors - 1), :) = z_traded;

% Compute z of services in 1972 in all countries
% z(:, i_services, :) = ...
%     compute_z_services(va, psi, p_sectoral, parameters);

z_services = zeros(n_countries, n_years);

for n = 1:n_countries
    for t = 1:n_years
        z_services(n, t) = ...
            xi^theta * ...
            B_replace(i_services, t)^theta * ...
            (va(n, i_services, t) / psi(n, i_services, t))^...
                                            (theta * beta_replace(i_services)) * ...
            prod(p_sectoral(n, :, t).^(gammas_replace(:, i_services, t)'))^theta * ...
            p_sectoral(n, i_services, t)^(-theta);
        
%         factor1(n, i_services, t) = theta * log(xi * B(i_services, t));
%         factor2(n, i_services, t) = theta * beta(i_services) * (log(va(n, i_services, t)) - log(psi(n, i_services, t)));
%         factor3(n, i_services, t) = 0;
%         factor4(n, i_services, t) = - theta * log(p_sectoral(n, i_services, t));
%         factor5(n, i_services, t) = theta * log(p_sectoral(n, :, t)) * gammas(:, i_services, t);
        
%         factor1(n, i_services, t) + factor2(n, i_services, t) + factor3(n, i_services, t) + factor4(n, i_services, t) + factor5(n, i_services, t)
        
    end % t
end % n

% for n = 1:n_countries
%     for t = 1:n_years
%         z_services(n, t) = ...
%             xi^theta * ...
%             B(i_services, t)^theta * ...
%             (va(n, i_services, t) / psi(n, i_services, t))^...
%                                             (theta * beta(i_services,t)) * ...
%             prod(p_sectoral(n, :, t).^(Gammas(:, i_services, t)'))^theta * ...
%             p_sectoral(n, i_services, t)^(-theta);
%         
% %         factor1(n, i_services, t) = theta * log(xi * B(i_services, t));
% %         factor2(n, i_services, t) = theta * beta(i_services) * (log(va(n, i_services, t)) - log(psi(n, i_services, t)));
% %         factor3(n, i_services, t) = 0;
% %         factor4(n, i_services, t) = - theta * log(p_sectoral(n, i_services, t));
% %         factor5(n, i_services, t) = theta * log(p_sectoral(n, :, t)) * gammas(:, i_services, t);
%         
% %         factor1(n, i_services, t) + factor2(n, i_services, t) + factor3(n, i_services, t) + factor4(n, i_services, t) + factor5(n, i_services, t)
%         
%     end % t
% end % n

z(:, i_services, :) = z_services;

%%
% f1 = factor2;
% f2 = f1 + factor3;
% f3 = f2 + factor4;
% f4 = f3 + factor5;
% f5 = f4 + factor1;
% 
% [~, f1] = detrend_series(f1, weights);
% [~, f2] = detrend_series(f2, weights);
% [~, f3] = detrend_series(f3, weights);
% [~, f4] = detrend_series(f4, weights);
% [~, f5] = detrend_series(f5, weights);
% 
% 
% n = 25;
% ff = [];
% ff(:, 1) = mean(corrcoef(squeeze(f1(n, :, :))'))';
% ff(:, 2) = mean(corrcoef(squeeze(f2(n, :, :))'))';
% ff(:, 3) = mean(corrcoef(squeeze(f3(n, :, :))'))';
% ff(:, 4) = mean(corrcoef(squeeze(f4(n, :, :))'))';
% ff(:, 5) = mean(corrcoef(squeeze(f5(n, :, :))'))';

% f = [];
% f(:, 1) = mean(corrcoef(squeeze(f1(n, :, :))'))';
% f(:, 2) = mean(corrcoef(squeeze(f2(n, :, :))'))';
% f(:, 3) = mean(corrcoef(squeeze(f3(n, :, :))'))';
% f(:, 4) = mean(corrcoef(squeeze(f4(n, :, :))'))';
% f(:, 5) = mean(corrcoef(squeeze(f5(n, :, :))'))';
% 
% 
% bb = {'va', 'trade', 'output p', 'input p', 'B'};
% 
% for aa = 1:5
%     figure(aa)
%     bar([f(:, aa), ff(:, aa)])
%     title([bb(aa), ' - Average correlation coefficient for all sectors - US' ])
%     xlabel('Sector')
%     legend('old', 'new', 'Location', 'SW')
% end
% 
% bar(ff)
% title('Average correlation coefficient for all sectors - US')
% xlabel('Sector')
% legend('va', 'trade', 'output p', 'input p', 'B', 'Location', 'SW')

% 
% bar(mean(corrcoef(squeeze(z(n, :, :))')))


% [~, ff] = detrend_series(log(z), weights);
% 
% zz = squeeze(ff(25, :, :));
% bar(mean(corrcoef(zz')))
% title('Average correlation coefficient for all sectors - US')
% xlabel('Sector')


%%
% zz = factor1 + factor2 + factor3 + factor4 + factor5;
% 
% [f1_t, f1_c] = detrend_series(factor1, weights);
% [f2_t, f2_c] = detrend_series(factor2, weights);
% [f3_t, f3_c] = detrend_series(factor3, weights);
% [f4_t, f4_c] = detrend_series(factor4, weights);
% [f5_t, f5_c] = detrend_series(factor5, weights);
% 
% [zz_t, zz_c] = detrend_series(zz, weights);
% 
% plot([squeeze(factor3(1, 1, :)), squeeze(f3_t(1, 1, :))])
% 
% 
% n = 25;
% for j = 1:24
%     aa = cov([squeeze(zz_c(n, j, :)), squeeze(f1_c(n, j, :)), squeeze(f2_c(n, j, :)), squeeze(f3_c(n, j, :)), squeeze(f4_c(n, j, :)), squeeze(f5_c(n, j, :))]);
%     
%     cc(j, :) = aa(1, 3:6) / aa(1, 1);
%     
%     bb = corrcoef([squeeze(f2_c(n, j, :)), squeeze(f3_c(n, j, :)), squeeze(f4_c(n, j, :)), squeeze(f5_c(n, j, :))]);
%     
%     dd(:, :, j) = bb;
%     
% end % for j
% 
% 
% mean(cc)
% 
% ee = mean(dd, 3);
% ff = mean(dd(:, :, 1:23), 3);
% ee(isnan(ee)) = ff(isnan(ee))


%% Compute Equipped Labor - L - still Cobb-Douglas
% This is alpha_j * (37) summed over j and then applying (38).
full_power = permute(repmat(alpha_replace ./ repmat(beta_replace * theta, [1 n_years]),...
                            [1, 1, n_countries]),...
                     [3 1 2]);
L = prod(z.^full_power, 2);
L = squeeze(L);



%% Save computed variables and parameters


baseline = struct;
load([c.model_folder, 'specification.mat'], 'model');
baseline.scenario = model;

baseline.country_names = country_names;
baseline.B = B_replace;
% baseline.K = K;
baseline.alpha = alpha_replace;
baseline.beta = beta_replace;
baseline.kappa = kappa;
baseline.theta = theta;
baseline.xi = xi;
baseline.gammas = gammas_replace;


deflator = pwt .* repmat(p_base, [1, n_countries]);
% baseline.trade_balance = bsxfun(@rdivide, trade_balance, p_base');
% baseline.trade_balance = trade_balance ./ deflator';
baseline.trade_balance = trade_balance;

baseline.L = L;
baseline.z = z;
baseline.psi = psi;

save([c.model_folder, 'alg_inputs.mat'], 'baseline');



%% Compute and save Real GDP and volatility in data

real_gdp_total = (va_total ./ deflator)';
real_gdp_sectoral = va ./ permute(repmat(deflator', [1 1 n_sectors]), [1 3 2]);

[~, data_cycle_total] = detrend_series(log(real_gdp_total), weights);
data_volatility_total = var(data_cycle_total, 0, 2);

save([c.model_folder, 'data_rgdp_and_volatility.mat'], 'country_names',...
     'data_volatility_total', 'deflator', 'pwt', 'p_base',...
     'va_total', 'va', 'p_sectoral_base', 'p_sectoral_data', 'p_sectoral', 'd');
end