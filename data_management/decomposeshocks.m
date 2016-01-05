function [Z_no_sector, Z_no_sector_residual] = decomposeshocks(Z, weights)

K = length(weights) - 1;

if K > 0
    % decompose log(Z) to trend and cycle
    [logZ_trend, logZ_cycle] = detrendseries(log(Z), weights);
    Z_hat = logZ_cycle;   
elseif K == 0;
    Z_hat = diff(log(Z), 1, 3); 
else
    error('Incorrect weights for filtering.')
end

% (Time) demean the growth rate series
Z_hat_mean = repmat(mean(Z_hat, 3), [1, 1, size(Z_hat, 3)]);
Z_tilde = Z_hat - Z_hat_mean;

% Sector-time effect: take average over countries
lambda = squeeze(mean(Z_tilde, 1));
save('results/lambdas.mat', 'lambda')
clear lambda
lambda_full = repmat(mean(Z_tilde, 1), [size(Z_tilde, 1), 1, 1]);

% Country-time effect: subtract sector-time specific effect and take average
% over sectors
% mu = squeeze(mean(Z_tilde - lambda_full, 2));
mu_full = repmat(mean(Z_tilde - lambda_full, 2), [1, size(Z_tilde, 2), 1]);

% Residual term
epsilon = Z_tilde - lambda_full - mu_full;

% calculate counterfactual productivities
Z_hat_no_sector = mu_full + epsilon + Z_hat_mean;
Z_hat_no_sector_residual = mu_full + Z_hat_mean;

if K > 0
    Z_no_sector = ...
        exp(logZ_trend + Z_hat_no_sector);
    Z_no_sector_residual = ...
        exp(logZ_trend + Z_hat_no_sector_residual);
else
    Z_no_sector = rollout(Z_hat_no_sector, Z);
    Z_no_sector_residual = ...
        rollout(Z_hat_no_sector_residual, Z);
end % if


% close all
% nn = 14;
% jj = 1;
% figure()
% subplot(3, 3, 1)
% plot(squeeze(Z(nn, jj, :)))
% title(['$Z_{' num2str(nn) ',t}^{' num2str(jj) '}$'], 'Interpreter', 'latex', 'FontSize', 15)
% %
% subplot(3, 3, 2)
% plot(squeeze(Z_hat(nn, jj, :)))
% title(['$\hat{Z}_{' num2str(nn) ',t}^{' num2str(jj) '}$ -- growth rate of Z'], 'Interpreter', 'latex', 'FontSize', 15)
% %
% subplot(3, 3, 3)
% plot(squeeze(Z_tilde(nn, jj, :)))
% title(['$\tilde{Z}_{' num2str(nn) ',t}^{' num2str(jj) '}$ -- demeaned growth rate of Z'], 'Interpreter', 'latex', 'FontSize', 15)
% %
% min_y = min([min(squeeze(lambda_full(nn, jj, :))), min(squeeze(mu_full(nn, jj, :))), min(squeeze(epsilon(nn, jj, :)))]);
% max_y = max([max(squeeze(lambda_full(nn, jj, :))), max(squeeze(mu_full(nn, jj, :))), max(squeeze(epsilon(nn, jj, :)))]);
% %
% subplot(3, 3, 4)
% plot(squeeze(lambda_full(nn, jj, :)))
% title(['$\lambda_{t}^{' num2str(jj) '}$ -- sectoral part'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 5)
% plot(squeeze(mu_full(nn, jj, :)))
% title(['$\mu_{' num2str(nn) ',t}$ -- country part'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 6)
% plot(squeeze(epsilon(nn, jj, :)))
% title(['$\epsilon_{' num2str(nn) ',t}^{' num2str(jj) '}$ -- residual part'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% min_y = min([min(log(squeeze(Z(nn, jj, :)))), min(log(squeeze(Z_no_sector(nn, jj, :)))), min(log(squeeze(Z_no_sector_residual(nn, jj, :)))) ]);
% max_y = max([max(log(squeeze(Z(nn, jj, :)))), max(log(squeeze(Z_no_sector(nn, jj, :)))), max(log(squeeze(Z_no_sector_residual(nn, jj, :)))) ]);
% %
% subplot(3, 3, 7)
% plot(log(squeeze(Z(nn, jj, :))))
% title(['$lnZ_{' num2str(nn) ',t}^{' num2str(jj) '}$'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 8)
% plot(log(squeeze(Z_no_sector(nn, jj, :))))
% title(['$lnZ_{' num2str(nn) ',t}^{' num2str(jj) '}$ no sectoral shocks'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 9)
% plot(log(squeeze(Z_no_sector_residual(nn, jj, :))))
% title(['$lnZ_{' num2str(nn) ',t}^{' num2str(jj) '}$ no sectoral and residual shocks'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 

end