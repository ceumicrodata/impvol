function [z_no_sector, z_no_sector_residual, z_mu_lambda] = decompose_shocks(z, weights, sectoral_weights)

assert(sum(sectoral_weights)==1);

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
dlmwrite('zhat.csv', z_hat, ',');
dlmwrite('ztilde.csv', z_tilde, ',');

% Sector-time effect: take average over countries
lambda = squeeze(mean(z_tilde, 1));
save('results/lambdas.mat', 'lambda')
clear lambda
lambda_full = repmat(mean(z_tilde, 1), [size(z_tilde, 1), 1, 1]);
dlmwrite('lambda.csv', lambda_full, ',');


% Country-time effect: subtract sector-time specific effect and take average
% over sectors
% weighted by sectoral weights
sw_full = repmat(sectoral_weights, [size(z_tilde, 1), 1, size(z_tilde, 3)]);
mu_full = repmat(sum(sw_full .* (z_tilde - lambda_full), 2), [1, size(z_tilde, 2), 1]);
dlmwrite('mu.csv', mu_full, ',');

% Residual term
epsilon = z_tilde - lambda_full - mu_full;
dlmwrite('epsilon.csv', epsilon, ',');

% calculate counterfactual productivities
z_hat_no_sector = mu_full + epsilon + z_hat_mean;
z_hat_no_sector_residual = mu_full + z_hat_mean;
z_mu_lambda = mu_full + lambda_full + z_hat_mean;

% verify that decomposition is exact
assert_all(...
    abs(log(z) - ...
        (logz_trend + z_hat_mean + mu_full + lambda_full + epsilon)) ...
    < 1e-5 ...
    )

if K > 0
    z_no_sector = ...
        exp(logz_trend + z_hat_no_sector);
    z_no_sector_residual = ...
        exp(logz_trend + z_hat_no_sector_residual);
    z_mu_lambda = ...
        exp(logz_trend + z_mu_lambda);
else
    z_no_sector = rollout(z_hat_no_sector, z);
    z_no_sector_residual = ...
        rollout(z_hat_no_sector_residual, z);
    z_mu_lambda = rollout(z_mu_lambda, z);
end % if

end