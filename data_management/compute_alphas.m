function alpha = compute_alphas(va, beta, gammas, weights)

[~, n_sectors, n_years] = size(va);

assert(n_sectors == size(gammas, 1), 'inconsistent number of sectors')

alpha = zeros(n_sectors, n_years);

for t = 1:n_years
    va_t = sum(va(:, :, t), 1)';
    alpha(:, t) = (eye(n_sectors) - gammas) * diag(1 ./ beta) * va_t / sum(va_t);
end % for t

% clip negativa alphas to zero
alpha(alpha < 0) = 0;

% smooth the series
[alpha, ~] = detrend_series(alpha, weights);

% normalize so that alphas sum to one in every period
alpha = bsxfun(@rdivide, alpha, sum(alpha, 1));

