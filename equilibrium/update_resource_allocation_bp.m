function L_share_njt = update_resource_allocation_bp(value_added_share, sectoral_wage_gap)

global c
weights = c.filter_weights;

[value_added_share_trend, ~] = ...
    detrend_series(value_added_share, weights);

L_share_star = bsxfun(@rdivide, value_added_share_trend, sum(value_added_share_trend, 2));

% if rho<infty, value added share should differ from ex ante optimal
if c.lac == 0
    L_share_njt = L_share_star;
else
    L_share_njt = L_share_star + (1 / c.lac) * sectoral_wage_gap; 
end
    
% adjust so that expectations are non-negative and sum up to one in each country 
L_share_njt = max(L_share_njt, c.numerical_zero);
L_share_njt = bsxfun(@rdivide, L_share_njt, sum(L_share_njt, 2));
end