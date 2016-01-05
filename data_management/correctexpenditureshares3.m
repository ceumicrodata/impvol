function d = correctexpenditureshares3(d, parameters)
% There are a number of manipulations to the shares...

global c

weights = c.filterWeights;

numericalZero = parameters.numericalZero;

[nCountries, ~, nSectors, nYears] = size(d);

% There is no data about domestic shares. We calculate them so that
% shares add up to 1.
for n = 1:nCountries
    d(n, n, :, :) = ones(1, 1, nSectors, nYears) - sum(d(n, :, :, :), 2);
end

domestic_shares = zeros(nCountries, nSectors - 1, nYears);

for t = 1:nYears
    for j = 1:(nSectors - 1)
        domestic_shares(:, j, t) = diag(d(:, :, j, t));
    end % j
end % t

% [domestic_trend, domestic_noise] = detrendseries(domestic_shares, weights);

% indd = domestic_shares < 0;
% t_c = squeeze(sum(indd, 2))';


% % cutoff = -0.126;
% cutoff = numericalZero;
% include_index = (domestic_trend(:) >= cutoff);
% exclude_index = logical(1 - include_index);
% 
% positive_index = (domestic_trend(:) > numericalZero);
% negative_index = logical(1 - positive_index);
% 
% 
% [d1, d2, d3] = size(domestic_trend);
% 
% % d1 = 2;
% % d2 = 3;
% % d3 = 4;
% country_D = repmat(eye(d1), [d2 * d3, 1]);
% sector_D = repmat(kron(eye(d2), ones(d1, 1)), [d3, 1]);
% year_D = kron(eye(d3), ones(d1 * d2, 1));
% country_D(:, 1) = [];
% sector_D(:, 1) = [];
% year_D(:, 1) = [];
% 
% country_sector_D = zeros(d1 * d2 * d3, (d1 - 1) * (d2 - 1));
% col = 0;
% for n = 1:(d1 - 1)
%     for j = 1:(d2 - 1)
%         col = col + 1;
%         country_sector_D(:, col) = country_D(:, n) .* sector_D(:, j);
%     end %j
% end % n
% 
% country_year_D = zeros(d1 * d2 * d3, (d1 - 1) * (d3 - 1));
% col = 0;
% for n = 1:(d1 - 1)
%     for t = 1:(d3 - 1)
%         col = col + 1;
%         country_year_D(:, col) = country_D(:, n) .* year_D(:, t);
%     end %j
% end % n
% 
% sector_year_D = zeros(d1 * d2 * d3, (d2 - 1) * (d3 - 1));
% col = 0;
% for j = 1:(d2 - 1)
%     for t = 1:(d3 - 1)
%         col = col + 1;
%         sector_year_D(:, col) = sector_D(:, j) .* year_D(:, t);
%     end %j
% end % n
% 
% c = ones(d1 * d2 * d3, 1);
% X_full = [c, country_D, sector_D, year_D, country_sector_D, country_year_D, sector_year_D];
% % rank(X_full)
% y_full = domestic_trend(:);
% 
% y = y_full(include_index);
% X = X_full(include_index, :);
% 
% XX_pinv = pinv(X' * X);
% next = XX_pinv * X';
% b_FE = next * y;
% y_hat = X * b_FE;
% y_hat_full = X_full * b_FE;
% % e = y - y_hat;
% 
% % y_imputed = y_full .* include_index + y_hat_full .* exclude_index;
% % y_imputed = max(min(y_imputed, 1), numericalZero);
% % domestic_trend_imputed = reshape(y_imputed, [nCountries, nSectors - 1, nYears]);
% domestic_trend_predicted = reshape(y_hat_full, [nCountries, nSectors - 1, nYears]);
% domestic_share_predicted = domestic_trend_predicted + domestic_noise ./ domestic_trend;

negative_indicator = squeeze(sum((domestic_shares < 0), 3));

domestic_shares_corrected = domestic_shares;

for nn = 1:nCountries
    for jj = 1:(nSectors - 1)
        if negative_indicator(nn, jj) > 0
%             pic_id = figure()
%             s2 = squeeze(domestic_trend(nn, jj, :));
            s1 = squeeze(domestic_shares(nn, jj, :));
            % s3 = squeeze(domestic_trend_predicted(nn, jj, :));
            s4 = min(max(s1, numericalZero), 1);
            s5 = ((s1 / (max(s1) - min(s1)) - min(s1) / (max(s1) - min(s1)))) * (max(max(s1), numericalZero) - max(min(s1), 0)) + max(min(s1), numericalZero);
            lam = 0.25;
            s6 = (1 - lam) * s4 + lam * s5;
            domestic_shares_corrected(nn, jj, :) = min(max(s6, numericalZero), 1);
%             plot(s1, '-k', 'LineWidth',2)
%             hold on
%             plot(s6, '-r', 'LineWidth',2)
%             grid on
%             set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
%             export_fig('results\\corrections.pdf', '-append')
%             close()
        end %if
    end % j
end % n
        
     
[domestic_shares_corrected, ~] = detrendseries(domestic_shares_corrected, weights);
domestic_shares_corrected = min(max(domestic_shares_corrected, numericalZero), 1);

for t = 1:nYears
    for j = 1:(nSectors - 1)
        d(:, :, j, t) = d(:, :, j, t) - diag(diag(d(:, :, j, t)));
        total_import_share = sum(d(:, :, j, t), 2) + numericalZero;
        total_import_share_imputed = 1 - domestic_shares_corrected(:, j, t); 
        d(:, :, j, t) = ...
            d(:, :, j, t) .* ...
            repmat(total_import_share_imputed ./ ...
                   total_import_share, [1, nCountries]);        
        d(:, :, j, t) = ...
            d(:, :, j, t) + diag(domestic_shares_corrected(:, j, t));
    end % j
end % t

% We will take logs and divide by shares, so we set a numerical zero
d(:, :, 1:(nSectors - 1), :) = ...
    max(d(:, :, 1:(nSectors - 1), :), numericalZero);
end