function [w_nj_new, P_nj_new, price_iterations, P_n_new, d] = wage_update(w_nj, L_nj, z_nj, P_nj, t)
% This function is one iteration of the wage loop.
% First it calculates prices corresponding to current wages,
% then calculates new wages.

global alpha beta theta kappa gammas S

[N, J] = size(L_nj);

% Compute D to be used in the price equations. See the definition of D in the notes.
% D depends only on current sector specific wages and aggregate labor

D = bsxfun(@times, permute((bsxfun(@power, bsxfun(@times, sum(L_nj, 2), w_nj)', -beta) .^ theta) .* z_nj', [3 2 1]), kappa(:, :, :, t).^theta);

% D2 = zeros(N, N, J);
% for j = 1:J
%     D2(:, :, j) = ...
%         repmat(z_nj(:, j)', N, 1) .* ...
%         (repmat((L_n .* w_nj(:, j))'.^(- beta(j)), N, 1) .* ...
%         kappa(:, :, j, t)).^theta;
% end
% max(abs((D(:) - D2(:))))

% Get sectoral prices. Note that this step depends on current wages through D.
[P_nj_new, price_iterations] = get_prices(P_nj, D, t);

% Calculate aggregate prices
P_n_new = prod(bsxfun(@power, bsxfun(@rdivide, P_nj_new, alpha(:, t)'), alpha(:, t)'), 2);

% preallocate and compute matrix for d values to be used in the wage equations
d = zeros(N, N, J); 
for j = 1:J
    d(:, :, j) = D(:, :, j) .* repmat(prod(bsxfun(@power, P_nj_new.^(-theta), gammas(:, j, t)'), 2)', N, 1);
    d(:, :, j) = d(:, :, j) ./ repmat(sum(d(:, :, j), 2), 1, N);
end

dd = reshape(permute(d, [3, 2, 1]), [N * J, N]);
R_jn = bsxfun(@rdivide, w_nj .* L_nj, beta')';
E_jn = (gammas(:, :, t) + alpha(:, t) * beta') * R_jn; % - alpha(:, t) * S(:, t)';
R_jn_new = reshape(sum(repmat(E_jn, [N, 1]) .* dd, 2), [J, N]);

% R_dif = 10;
% R_tol = 1e-4;
% while R_dif > R_tol
%     E_jn = (gammas + alpha(:, t) * beta') * R_jn - alpha(:, t) * S(:, t)' * 1000; % convert S to millions
%     R_jn_new = reshape(sum(repmat(E_jn, [N, 1]) .* dd, 2), [J, N]);
%     R_dif = max(abs(R_jn_new(:) - R_jn(:)));
%     R_jn = R_jn_new;
% end

w_nj_new = bsxfun(@times, R_jn_new' ./ L_nj, beta');

end
