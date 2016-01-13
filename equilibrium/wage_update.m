function [w_nj_new, P_n, price_iterations, price_lambda, P_nj, d] = wage_update(w_nj, L_nj, L_n, z_nj, P_n, t)
% This function is one iteration of the wage loop.
% First it calculates prices corresponding to current wages,
% then calculates new wages.

global alpha beta theta kappa

[N, J] = size(L_nj);

% Compute D to be used in the price equations. See the definition of D in the notes.
% D depends only on current sector specific wages and aggregate labor
D = zeros(N, N, J);
for j = 1:J
    D(:, :, j) = ...
        repmat(z_nj(:, j)', N, 1) .* ...
        (repmat((L_n .* w_nj(:, j))'.^(- beta(j)), N, 1) .* ...
        kappa(:, :, j, t)).^theta;
end

% D = bsxfun(@times, permute((bsxfun(@power, bsxfun(@times, sum(L_nj, 2), w_nj)', -beta) .^ theta) .* z_nj', [3 2 1]), kappa(:, :, :, t).^theta);

% Get aggregate prices. Note that this depends on current wages through D!
[Ptheta_n, price_iterations, price_lambda] = get_prices(D, P_n.^theta, t);
P_n = Ptheta_n .^ (1 / theta);

% Calculate sectoral prices
P_nj = recover_sectoral_prices(D, P_n);

% preallocate and compute matrix for d values to be used in the wage equations
d = zeros(N, N, J); 
for j = 1:J
    d(:, :, j) = D(:, :, j) .* repmat(Ptheta_n'.^(beta(j) - 1), N, 1);
    d(:, :, j) = d(:, :, j) ./ repmat(sum(d(:, :, j), 2), 1, N);
end

% AA = bsxfun(@times, D , permute(bsxfun(@power, repmat(Ptheta_n', [J 1]), beta - 1), [3 2 1]));
% d = bsxfun(@rdivide, AA, sum(AA, 2));

% compute A, the matrix in the wage equations
A1 = [];
A2 = [];
for j = 1:J
    A1 = vertcat(A1, alpha(j, t) * beta(j) * repmat(d(:, :, j)', 1, J));
    A2 = blkdiag(A2, (1 - beta(j)) * d(:, :, j)');
end
A = A1 + A2;

% Do one step in the iterative procedure
wL_nj = w_nj(:) .* L_nj(:);
wL_nj_new = A * wL_nj;
    
% calculate wages from value added by dividing by labor hours
w_nj_new = reshape(wL_nj_new, N, J) ./ L_nj;
end
