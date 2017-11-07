function compute_volatilities

global c

load([c.model_folder, 'equilibrium.mat'])

% sector, time, country
log_sectoral_GDP = permute(log(equilibrium.w_njt) + log(equilibrium.L_njt), [2 3 1]);

realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
volatilities = var(cycle, 0, 2);

country_names = importdata([c.data_folder_original, 'country_name.txt']);
% % China removed
% country_names(strcmp('China', country_names)) = [];
% %%%%%%%%%%%%%%%
fid = fopen([c.model_folder, 'volatilities.csv'], 'w');
fprintf(fid, '%25s, %12s\n',  'country', 'volatility');

for ii = 1:(length(volatilities) - 1)
    fprintf(fid, '%25s, %12.10f\n', country_names{ii}, volatilities(ii));
end    
fprintf(fid, '%25s, %12.10f', country_names{end}, volatilities(end));
fclose(fid);


