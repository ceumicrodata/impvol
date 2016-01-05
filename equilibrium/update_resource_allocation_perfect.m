function [L_njt_new] = update_resource_allocation_perfect(log_value_added_share, L_nt)


[~, J, ~] = size(log_value_added_share);

E_value_added_share = exp(log_value_added_share);
               
% multiply by total labor to get sectoral employment
L_nt_full = permute(repmat(L_nt, [1 1 J]), [1 3 2]);

L_njt_new = L_nt_full .* E_value_added_share;

% nn = 2;
% jj = 1;

% pic_id = figure();
% 
% plot(squeeze(exp(log_value_added_share(nn, jj, :))))
% ylim([0.005, 0.02])  
% print(pic_id, 'value_added_shares.pdf', '-dpsc', '-append')
% close(pic_id)


end