function [L_njt_new] = update_resource_allocation_bp(log_value_added_share, L_nt)

global c
weights = c.filter_weights;

[~, J, ~] = size(log_value_added_share);

[log_value_added_share_trend, ~] = ...
    detrend_series(log_value_added_share, weights);
   
% calculate expected wagebill shares
E_value_added_share = exp(log_value_added_share_trend);

% adjust so that expectations sum up to one in each country
E_value_added_share = E_value_added_share ./ ...
                   repmat(sum(E_value_added_share, 2), [1 J 1]);
               
% multiply by total labor to get sectoral employment
L_nt_full = permute(repmat(L_nt, [1 1 J]), [1 3 2]);

L_njt_new = L_nt_full .* E_value_added_share;
end