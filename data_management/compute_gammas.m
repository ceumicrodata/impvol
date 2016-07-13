function gammas = compute_gammas(io_values, total_output, output_shares, intermediate_input_shares)

years = 13;
sectors = 34;

io_values = permute(reshape(io_values, [sectors, sectors, years]), [2, 1, 3]);
total_output = reshape(total_output, [sectors, years]);

% add up agriculture and mining
io_values(1, :, :) = sum(io_values(1:2, :, :));
io_values(2, :, :) = [];

io_values(:, 1, :) = sum(io_values(:, 1:2, :), 2);
io_values(:, 2, :) = [];

total_output(1, :) = sum(total_output(1:2, :));
total_output(2, :) = [];

% add up services
io_values(18, :, :) = sum(io_values(18:end, :, :));
io_values(19:end, :, :) = [];

io_values(:, 18, :) = sum(io_values(:, 18:end, :), 2);
io_values(:, 19:end, :) = [];

total_output(18, :) = sum(total_output(18:end, :));
total_output(19:end, :) = [];



% split rows
% duplicate the rows that will be split
dupl_index = [1, 2, 2, 3, 3, 3, 4, 5, 5, 6:12, 13, 13, 13, 14:18]';
io_values_dupl = io_values(dupl_index, :, :);

% set the shares used in the split to 1
output_shares_full = ones(size(io_values_dupl));

% set the shares of the affected sectors to the appropriate values
split_index = [2:6, 8:9, 17:19]';
output_shares_full(split_index, :, :) = ...
    permute(repmat(output_shares, [1, 1, size(io_values_dupl, 2)]), [2, 3, 1]);

io_values_new = io_values_dupl .* output_shares_full;


% split columns
% duplicate the columns that will be split
io_values_dupl = io_values_new(:, dupl_index, :);

% output_shares_full = ones(size(io_values_dupl));
% output_shares_full(:, split_index, :) = ...
%     permute(repmat(output_shares, [1, 1, size(io_values_dupl, 1)]), [3, 2, 1]);
intermediate_input_shares_full = ones(size(io_values_dupl));
intermediate_input_shares_full(:, split_index, :) = ...
    permute(repmat(intermediate_input_shares, [1, 1, size(io_values_dupl, 1)]), [3, 2, 1]);

io_values_new = io_values_dupl .* intermediate_input_shares_full;
% io_values_new = io_values_dupl .* output_shares_full;



% split total output
total_output = total_output(dupl_index, :);
output_shares_full = ones(size(total_output));
output_shares_full(split_index, :) = output_shares';
total_output = total_output .* output_shares_full;



% correct the order of sectors
io_values_new([18, 19, 20], :, :) = io_values_new([20, 18, 19], :, :);
io_values_new(:, [18, 19, 20], :) = io_values_new(:, [20, 18, 19], :);
total_output([18, 19, 20], :) = total_output([20, 18, 19], :); 



% get gammas
gammas = bsxfun(@rdivide, io_values_new, permute(total_output, [3, 1, 2]));
gammas = mean(gammas, 3);