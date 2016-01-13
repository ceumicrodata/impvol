function [z_no_sector, z_no_sector_residual] = decompose_shocks(z, weights)

K = length(weights) - 1;

if K > 0
    % decompose log(z) to trend and cycle
    [logz_trend, logz_cycle] = detrend_series(log(z), weights);
    z_hat = logz_cycle;   
elseif K == 0;
    z_hat = diff(log(z), 1, 3); 
else
    error('Incorrect weights for filtering.')
end

% (Time) demean the growth rate series
z_hat_mean = repmat(mean(z_hat, 3), [1, 1, size(z_hat, 3)]);
z_tilde = z_hat - z_hat_mean;

% Sector-time effect: take average over countries
lambda = squeeze(mean(z_tilde, 1));
save('results/lambdas.mat', 'lambda')
clear lambda
lambda_full = repmat(mean(z_tilde, 1), [size(z_tilde, 1), 1, 1]);

% Country-time effect: subtract sector-time specific effect and take average
% over sectors
% mu = squeeze(mean(z_tilde - lambda_full, 2));
mu_full = repmat(mean(z_tilde - lambda_full, 2), [1, size(z_tilde, 2), 1]);

% Residual term
epsilon = z_tilde - lambda_full - mu_full;

% calculate counterfactual productivities
z_hat_no_sector = mu_full + epsilon + z_hat_mean;
z_hat_no_sector_residual = mu_full + z_hat_mean;

if K > 0
    z_no_sector = ...
        exp(logz_trend + z_hat_no_sector);
    z_no_sector_residual = ...
        exp(logz_trend + z_hat_no_sector_residual);
else
    z_no_sector = rollout(z_hat_no_sector, z);
    z_no_sector_residual = ...
        rollout(z_hat_no_sector_residual, z);
end % if


% close all
% nn = 14;
% jj = 1;
% figure()
% subplot(3, 3, 1)
% plot(squeeze(z(nn, jj, :)))
% title(['$z_{' num2str(nn) ',t}^{' num2str(jj) '}$'], 'Interpreter', 'latex', 'FontSize', 15)
% %
% subplot(3, 3, 2)
% plot(squeeze(z_hat(nn, jj, :)))
% title(['$\hat{z}_{' num2str(nn) ',t}^{' num2str(jj) '}$ -- growth rate of z'], 'Interpreter', 'latex', 'FontSize', 15)
% %
% subplot(3, 3, 3)
% plot(squeeze(z_tilde(nn, jj, :)))
% title(['$\tilde{z}_{' num2str(nn) ',t}^{' num2str(jj) '}$ -- demeaned growth rate of z'], 'Interpreter', 'latex', 'FontSize', 15)
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
% min_y = min([min(log(squeeze(z(nn, jj, :)))), min(log(squeeze(z_no_sector(nn, jj, :)))), min(log(squeeze(z_no_sector_residual(nn, jj, :)))) ]);
% max_y = max([max(log(squeeze(z(nn, jj, :)))), max(log(squeeze(z_no_sector(nn, jj, :)))), max(log(squeeze(z_no_sector_residual(nn, jj, :)))) ]);
% %
% subplot(3, 3, 7)
% plot(log(squeeze(z(nn, jj, :))))
% title(['$lnz_{' num2str(nn) ',t}^{' num2str(jj) '}$'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 8)
% plot(log(squeeze(z_no_sector(nn, jj, :))))
% title(['$lnz_{' num2str(nn) ',t}^{' num2str(jj) '}$ no sectoral shocks'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 
% %
% subplot(3, 3, 9)
% plot(log(squeeze(z_no_sector_residual(nn, jj, :))))
% title(['$lnz_{' num2str(nn) ',t}^{' num2str(jj) '}$ no sectoral and residual shocks'], 'Interpreter', 'latex', 'FontSize', 15)
% ylim([min_y max_y]) 

end