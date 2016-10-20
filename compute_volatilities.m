function compute_volatilities

global c

load([c.model_folder, 'equilibrium.mat'])

realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
volatilities = var(cycle, 0, 2);

country_names = importdata([c.data_folder_original, 'country_name.txt']);
fid = fopen([c.model_folder, 'volatilities.csv'], 'w');
fprintf(fid, '%25s, %12s\n',  'country', 'volatility');

for ii = 1:(length(volatilities) - 1)
    fprintf(fid, '%25s, %12.10f\n', country_names{ii}, volatilities(ii));
end    
fprintf(fid, '%25s, %12.10f', country_names{end}, volatilities(end));
fclose(fid);


