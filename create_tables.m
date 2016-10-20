function create_tables(model_list_file, table_name)

model_list = importdata(model_list_file);

% load volatilities from model folders
volatilities = [];
for ii = 1:4
    model = model_list{ii};
    model_vol = importdata(['models/', model, '/volatilities.csv']);
    volatilities(:, ii) = model_vol.data;
end

% compute new columns in table
new_cols = [];
new_cols(:, 1) = 100 * (volatilities(:, 1) ./ volatilities(:, 3) - 1);
new_cols(:, 3) = 100 * (volatilities(:, 2) - volatilities(:, 4)) ./ volatilities(:, 3);
new_cols(:, 2) = new_cols(:, 1) - new_cols(:, 3);

% load data volatilities
model = model_list{1};
load(['models/', model, '/data_rgdp_and_volatility.mat'], 'data_volatility_total');

fid = fopen(table_name, 'w');
for ii = 1:size(volatilities, 1)
    fprintf(fid, '%25s, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.4f, %12.4f, %12.4f\n',...
        model_vol.textdata{ii + 1, 1}, data_volatility_total(ii), volatilities(ii, :), new_cols(ii, :));
end
fclose(fid);




