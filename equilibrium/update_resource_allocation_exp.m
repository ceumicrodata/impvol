function [L_njt_new] = update_resource_allocation_exp(log_value_added_share, ...
                                                      L_nt)
                               
[N, J, T] = size(log_value_added_share);

y = reshape(log_value_added_share, [N * J, T])';

% generate regressors: a constant and a time trend
X = [ones(T, 1), (1:T)'];
% calculate the projection matrix
P_X = (X'*X) \ X';

% calculate all coefficients
delta = P_X * y;
% delta(:, 1:25)

% calculate the fitted values
y_hat = X * delta;
    
% calculate expected wagebill shares
E_value_added_share = reshape(exp(y_hat'), [N J T]);
    
% adjust so that expectations sum up to one in each country
E_value_added_share = E_value_added_share ./ ...
                   repmat(sum(E_value_added_share, 2), [1 J 1]);
    
% multiply by total labor to get sectoral employment
L_nt_full = permute(repmat(L_nt, [1 1 J]), [1 3 2]);

L_njt_new = L_nt_full .* E_value_added_share;

nn = 13;
jj = 21;

pic_id = figure();

plot([squeeze(exp(log_value_added_share(nn, jj, :))),...
      squeeze(E_value_added_share(nn, jj, :))])
ylim([0, 0.2])  
print(pic_id, 'value_added_shares.pdf', '-dpsc', '-append')
close(pic_id)

% save the plot of fitted values for the same sector and country in every iteration
%pic_id = figure(outer_iteration);
%figure()
%plot([y(:, 600), X * delta(:, 600)])
%plot([X * delta(:, sector_to_plot, country_to_plot), y(sector_to_plot, :, country_to_plot)'])
%titlestring1 = 'Wagebill time series and the fitted regression line';
%titlestring2 = strjoin({'Country:', strjoin(names(country_to_plot)) ,...
%                        'Sector:', num2str(sector_to_plot),...
%                        'Iteration:' , num2str(outer_iteration)});
%title({titlestring1; titlestring2})
%xlabel('Time period')
%ylabel('Log Wagebill Share')
%print(pic_id, 'regression_figures_cf2.pdf', '-dpsc', '-append')
%close(pic_id)

end