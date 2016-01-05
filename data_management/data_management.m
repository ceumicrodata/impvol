function data_management
%% Data management

global c

FE = c.FE;
FE_prices = c.FE_prices;
alphaType = c.alphaType;
iBase = c.iBase;
havePrices = c.havePrices;

theta = c.theta;
eta = c.eta;
numericalZero = c.numericalZero;
weights = c.filterWeights;

inputFolder = c.dataFolderOriginal;
outputFolder = c.dataFolderAlgorithmInput;

parameters.theta = theta;
parameters.numericalZero = numericalZero;
parameters.iBase = iBase;



%% Import files from data folder.
countryNames = importdata([inputFolder, 'country_name.txt']);
betaPanel = dlmread([inputFolder, 'beta_panel.txt'], '\t', 1, 2);

pwt = dlmread([inputFolder,...
    'aggregate_price_relative_to_US.csv'], ',', 0, 2);

va = dlmread([inputFolder, 'sectoral_value_added.txt'], '\t', 1, 2);
importShares = dlmread([inputFolder, 'import_share.txt'], '\t');


%% Get key constants of data dimensions
nCountries = length(countryNames);
nSectors = size(betaPanel, 2);
nYears = size(betaPanel, 1) / nCountries;
iServices = nSectors;

parameters.nCountries = nCountries;
parameters.nSectors = nSectors;
parameters.nYears = nYears;
parameters.iServices = iServices;



%% Compute alphas
% Reshape value added to a more convenient 3 dimensional form with
% dimensions (country, sector, year).
va = reshape(reshape(va, nYears, nCountries * nSectors)',...
    nCountries, nSectors, nYears);

% To get sectoral betas take average across years and countries.
beta = mean(betaPanel)';

if alphaType == 0
    % Constant alpha
    va_sum_over_n_t = squeeze(sum(sum(va, 1), 3));
    sectoral_va_share = va_sum_over_n_t / sum(va_sum_over_n_t);
    alpha = (sectoral_va_share' ./ beta) / sum(sectoral_va_share' ./ beta);
    alpha = repmat(alpha, [1 nYears]);
else
    % Time varying alpha
    va_sum_over_n = squeeze(sum(va, 1));
    sectoral_va_share = va_sum_over_n ./ ...
                        repmat(sum(va_sum_over_n, 1), [nSectors 1]);
    sectoral_va_share_over_beta = sectoral_va_share ./ ...
                                  repmat(beta, [1 nYears]);
    
    alpha = sectoral_va_share_over_beta ./ ...
            repmat(sum(sectoral_va_share_over_beta, 1), [nSectors 1]);
      
    % Smooth the time varying alphas using a band-pass filter
    if alphaType == 2
        [alphaTrend, ~] = detrendseries(alpha, weights);
        alpha = alphaTrend ./ repmat(sum(alphaTrend, 1), [nSectors 1]);
    end % if
end % if



%% Import sectoral Prices
pSectoralData = csvread([inputFolder, 'sectoral_price_index.csv']);
ii = sum(havePrices(1:iBase));
pSectoralBase = pSectoralData((ii - 1) * nYears + 1 : ii * nYears, :);



%% Process imported matrices
% Normalize sectoral price index in base country 
pSectoralBase = pSectoralBase ./ repmat(pSectoralBase(1, :), [nYears 1]);

% Compute aggregate price index
pBase = prod((pSectoralBase' ./ alpha) .^ alpha)';

% Recompute PWT. Needed if US is not chosen as the base country
pwt = reshape(pwt, nYears, nCountries);
pwt = pwt ./ repmat(pwt(:, iBase), [1 nCountries]);

% Convert import shares to expenditure shares based on source country.
d = zeros(nCountries, nCountries, nSectors, nYears);
nLinesInImportShares = length(importShares);

firstYear = importShares(1, 1);
for iLine = 1:nLinesInImportShares
    year = importShares(iLine, 1) + 1 - firstYear;
    importer = importShares(iLine, 2);
    exporter = importShares(iLine, 3);
    
    % Clip negative import shares from below at 0.
    d(importer, exporter, 1:(nSectors - 1), year) = ...
        max(importShares(iLine, 4:26), 0);
    
    % Set services shares to 0 as they are not traded.
    d(importer, exporter, iServices, year) = 0;
end

% There are a number of manipulations to the shares...
% d1 = correctexpenditureshares(d, parameters);
% d2 = correctexpenditureshares2(d, parameters);
d = correctexpenditureshares(d, parameters);


% nn = 3;
% NN = 3;
% jj = 4;
% tt = 1:36;
% s1 = squeeze(d1(nn, NN, jj, tt));
% s2 = squeeze(d2(nn, NN, jj, tt));
% s3 = squeeze(d3(nn, NN, jj, tt));
% figure()
% plot([s1, s2, s3], 'LineWidth', 2)


%% Calculate new matrices of constants.
xi = gamma((theta + 1 - eta) / theta);

B = beta.^(-beta) .* (1 - beta).^(beta - 1);
K = zeros(nYears, 1);
for t = 1:nYears
    K(t) = prod((xi * B ./ alpha(:, t)).^(alpha(:, t)))^theta;
end % t

vaTotal = sum(va, 2);
vaTotal = squeeze(vaTotal)';

kappa = computetradecost(d, parameters);

% plot_trade_costs(kappa)

% collect parameters
parameters.alpha = alpha;
parameters.B = B;
parameters.beta = beta;
parameters.xi = xi;
parameters.kappa = kappa;



%% Get expectations of value added shares (psi)
vaShares = va ./ permute(repmat(vaTotal', [1 1 nSectors]), [1 3 2]);
[psi, ~] = detrendseries(vaShares, weights);



%% Calculate the log difference in variables
pSectoralBaseChange = diff(log(pSectoralBase));
pwtChange = diff(log(pwt));
pBaseChange = diff(log(pBase));

% vaTotalChange = diff(log(vaTotal));
vaChange = diff(log(va), 1, 3);
psiChange = diff(log(psi), 1, 3);

kappaToBaseChange = diff(log(squeeze(kappa(iBase, :, :, :))), 1, 3);
dToBaseChange = diff(log(squeeze(d(iBase, :, :, :))), 1, 3);



%%
if FE
    zeta = computeZeta(d, va, psi, pwt, pBase, parameters);
    
    %save('zeta.mat', 'zeta')
    
    tau_base_single = theta * log(pSectoralBase(:, 1:(nSectors - 1))');
    tau_base = permute(repmat(tau_base_single, [1 1 nCountries]), [3 1 2]);
    
    zeta_base = zeta(iBase, :, :, :);
    zeta_base = repmat(zeta_base, [nCountries 1 1 1]);
    
    tau = squeeze(mean(zeta - zeta_base, 2)) + tau_base;
    tau_full = permute(repmat(tau, [1 1 1 nCountries]), [1 4 2 3]);
    
    lambda = squeeze(mean(zeta - tau_full, 1));
    
    Z_traded = exp(lambda);
    Z = zeros(nCountries, nSectors, nYears);
    Z(:, 1:(nSectors - 1), :) = Z_traded;
    
    pSectoral = zeros(nCountries, nSectors, nYears);
    pSectoral(:, 1:(nSectors - 1), :) = exp(tau / theta);
    
else
    % Compute Shocks - except for services
    % Compute the normalized aggregate price in the US in 1972
    pBase1972 = pBase(1);
    parameters.pBase1972 = pBase1972;
    
    noBaseIndex = [1:nCountries];
    noBaseIndex(iBase) = [];
    parameters.noBaseIndex = noBaseIndex;
    
    % Initialize Z
    Z = zeros(nCountries, nSectors, nYears);
    
    % Compute Z in 1972 for all other countries based on (34) and (35)
    % Note: Z of services is computed only for the base country
    Z(:, :, 1) = computeZ1972(d, va, psi, pwt, parameters);
    
    % Compute change of Z(US, j, t) based on (33) and (31)
    % Note: Change of Z for services is computed only for the base country
    changeZ = computechangeZ(dToBaseChange, vaChange, psiChange,...
        pBaseChange, pSectoralBaseChange, pwtChange,...
        kappaToBaseChange, parameters);
    
    % Compute Z(n, j, t) except for services from changes and init. values
    Z = rollout(changeZ, Z);
end
% 
% 
% 
%% Compute Shocks - for services in all countires other than the US
% Compute Phi(n, j, t), except for services
Phi = computePhi(Z, va, psi, pwt, pBase, d, parameters);
%
% Compute sectoral prices
if ~FE_prices
    pSectoral = zeros(nCountries, nSectors, nYears);
    pSectoral(:, 1:(nSectors - 1), :) = xi * (Phi).^(- 1 / theta);
end

pSectoral(:, iServices, :) = ...
    computepServices(pwt, pSectoral, pBase, parameters);
% 
% % Compute Z of services in 1972 in all countries
Z(:, iServices, :) = ...
    computeZServices(va, psi, pwt, pBase, pSectoral, parameters);

%save('sectoral_prices', 'pSectoral')

%% Compute Equipped Labor - L
% This is alpha_j * (37) summed over j and then applying (38).
full_power = permute(repmat(alpha ./ repmat(beta * theta, [1 nYears]),...
                            [1, 1, nCountries]),...
                     [3 1 2]);
L = prod(Z.^full_power, 2);
L = squeeze(L);



%% Save computed variables and parameters
baseline = struct;
baseline.scenario = 'Baseline';

baseline.countryNames = countryNames;
baseline.B = B;
baseline.K = K;
baseline.alpha = alpha;
baseline.beta = beta;
baseline.kappa = kappa;
baseline.theta = theta;
baseline.xi = xi;

baseline.L = L;
baseline.Z = Z;
baseline.psi = psi;

save([outputFolder, 'data_theta_', num2str(theta), '.mat'], 'baseline');



%% Compute and save Real GDP and volatility in data
% deflator = pwt .* repmat(pUs * pUs1972, [1, nCountries]);
deflator = pwt .* repmat(pBase, [1, nCountries]);

realGdpTotal = (vaTotal ./ deflator)';
realGdpSectoral = va ./ permute(repmat(deflator', [1 1 nSectors]), [1 3 2]);

[~, dataCycleTotal] = detrendseries(log(realGdpTotal), weights);
dataVolatilityTotal = var(dataCycleTotal, 0, 2);

outputFolder = c.resultsFolder;
save([outputFolder, 'data_rgdp_and_volatility.mat'], 'realGdpSectoral',...
     'realGdpTotal', 'dataVolatilityTotal', 'deflator', 'pwt', 'pBase',...
     'vaTotal', 'va', 'pSectoralBase', 'pSectoralData', 'd');

end