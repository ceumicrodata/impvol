function d = correct_expenditure_shares(d, parameters)
% There are a number of manipulations to the shares...

numerical_zero = parameters.numerical_zero;

[n_countries, ~, n_sectors, n_years] = size(d);

share_withing_imports = d./repmat(sum(d, 2), [1, n_countries, 1, 1]);
domestic_per_total_import = 1/repmat(sum(d, 2),[1, n_countries, 1, 1]) - 1;
domestic_per_total_import(domestic_per_total_import < numerical_zero) = numerical_zero;
d = share_withing_imports./ (1 + domestic_per_total_import);

% There is no data about domestic shares. We calculate them so that
% shares add up to 1.
for n = 1:n_countries
    d(n, n, :, :) = ones(1, 1, n_sectors, n_years) - sum(d(n, :, :, :), 2);
end

% (Potentially) We created some negative domestic shares that we have to
% correct for.

% We clip them from below to 0 and then renormalize the
% shares so that they add up to 1 for a given (destination) country, sector
% time period triplet.
%d(d < 0) = 0;
%d = d ./ repmat(sum(d, 2), [1, n_countries, 1, 1]);

% % In the first period we clip them from below to 0 and then renormalize the
% % shares so that they add up to 1 for a given (destination) country, sector
% % time period triplet.
% d(d(:, :, :, 1) < 0) = 0;
% d(:, :, :, 1) = ...
%     d(:, :, :, 1) ./ repmat(sum(d(:, :, :, 1), 2), [1, n_countries, 1]);
% 
% % In the later years we approximate the negative domestic share with its
% % previous value and then renormalize again.
% for t = 2:n_years
%     d_current_year = d(:, :, :, t);
%     d_previous_year = d(:, :, :, t - 1);
%     index_of_negative = (d_current_year < 0);
%     d_current_year(index_of_negative) = d_previous_year(index_of_negative);
%     d(:, :, :, t) = ...
%         d_current_year ./ repmat(sum(d_current_year, 2), [1, n_countries, 1]);
% end
% 
% We will take logs and divide by shares, so we set a numerical zero
d(:, :, 1:(n_sectors - 1), :) = ...
    max(d(:, :, 1:(n_sectors - 1), :), numerical_zero);
end
