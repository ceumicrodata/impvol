function create_table(table_name)

fid = fopen(['model_specifications/', table_name, '.csv']);
model_list = textscan(fid, '%s', 5, 'Delimiter', ',');
fclose(fid);

% load volatilities from model folders
volatilities = [];
for i = 1:4
    model = [table_name, '_', model_list{1}{i + 1}];
    model_vol = importdata(['models/', model, '/volatilities.csv']);
    volatilities(:, i) = model_vol.data;
end


% compute new columns in table
new_cols = [];
new_cols(:, 1) = 100 * (volatilities(:, 1) ./ volatilities(:, 3) - 1);
new_cols(:, 3) = 100 * (volatilities(:, 2) - volatilities(:, 4)) ./ volatilities(:, 3);
new_cols(:, 2) = new_cols(:, 1) - new_cols(:, 3);

% load data volatilities
model = model_list{1}{2};
load(['models/', table_name, '_', model, '/data_rgdp_and_volatility.mat'], 'data_volatility_total');

fid = fopen(['tables/', table_name, '.csv'], 'w');
fprintf(fid, 'country, data, baseline, nosectoral, kappa1972, kappa1972_nosectoral, trade_barriers, specialization, diversification\n');
for ii = 1:size(volatilities, 1)
    fprintf(fid, '%25s, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.4f, %12.4f, %12.4f\n',...
        model_vol.textdata{ii + 1, 1}, data_volatility_total(ii), volatilities(ii, :), new_cols(ii, :));
end
fclose(fid);