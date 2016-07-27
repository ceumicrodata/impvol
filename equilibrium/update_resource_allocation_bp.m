function E_value_added_share = update_resource_allocation_bp(log_value_added_share, sectoral_wage_gap)

global c
weights = c.filter_weights;

% if rho<infty, value added share should differ from ex ante optimal
% FIXME: what if this share becomes 0 or negative?
input_series = log(exp(log_value_added_share) - sectoral_wage_gap*c.labor_adjustment_cost);

[~, J, ~] = size(input_series);

[log_value_added_share_trend, ~] = ...
    detrend_series(input_series, weights);
   
% calculate expected wagebill shares
E_value_added_share = exp(log_value_added_share_trend);

% adjust so that expectations sum up to one in each country
E_value_added_share = (E_value_added_share ./ ...
                   repmat(sum(E_value_added_share, 2), [1 J 1])) ...
					+ sectoral_wage_gap*c.labor_adjustment_cost);
               
% multiply by total labor to get sectoral employment
% L_nt_full = permute(repmat(L_nt, [1 1 J]), [1 3 2]);

% L_njt_new = L_nt_full .* E_value_added_share;
end