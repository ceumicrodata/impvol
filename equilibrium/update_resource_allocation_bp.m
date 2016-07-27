function E_value_added_share = update_resource_allocation_bp(log_value_added_share, sectoral_wage_gap)

global c
weights = c.filter_weights;

% if rho<infty, value added share should differ from ex ante optimal
input_series = exp(log_value_added_share) - sectoral_wage_gap*c.labor_adjustment_cost;

[~, J, ~] = size(input_series);

[value_added_share_trend, ~] = ...
    detrend_series(input_series, weights);
   
% adjust so that expectations sum up to one in each country
E_value_added_share = (value_added_share_trend ./ ...
                   repmat(sum(value_added_share_trend, 2), [1 J 1])) ...
					+ sectoral_wage_gap*c.labor_adjustment_cost;
               
end