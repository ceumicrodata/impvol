function [w_nj_new, P_nj_new, price_iterations, P_n_new, d] = wage_update(w_nj, L_nj, z_nj, P_nj, t, va_to_fit, p_to_fit, B_gamma, B_beta, D_alpha, S_full, beta_full)
% This function is one iteration of the wage loop.
% First it calculates prices corresponding to current wages,
% then calculates new wages.

global c alpha beta theta kappa gammas 

i_base = c.i_base;

[N, J] = size(L_nj);

% Compute D to be used in the price equations. See the definition of D in the notes.
% D depends only on current sector specific wages and aggregate labor
D = bsxfun(@times, permute((bsxfun(@power, bsxfun(@times, sum(L_nj, 2), w_nj)', -beta) .^ theta) .* z_nj', [3 2 1]), kappa(:, :, :, t).^theta);

% Get sectoral prices. Note that this step depends on current wages through D.
[P_nj_new, price_iterations] = get_prices(P_nj, D, t);

% Calculate aggregate prices
P_n_new = prod(bsxfun(@power, bsxfun(@rdivide, P_nj_new, alpha(:, t)'), alpha(:, t)'), 2);

% xx = p_to_fit / P_n_new(i_base);
% 
% P_n_new = xx * P_n_new;
% P_nj_new = xx * P_nj_new;
% w_nj = xx * w_nj;
% 
% D = bsxfun(@times, permute((bsxfun(@power, bsxfun(@times, sum(L_nj, 2), w_nj)', -beta) .^ theta) .* z_nj', [3 2 1]), kappa(:, :, :, t).^theta);


% preallocate and compute matrix for d values to be used in the wage equations
d = zeros(N, N, J); 
for j = 1:J
    d(:, :, j) = D(:, :, j) .* repmat(prod(bsxfun(@power, P_nj_new, - theta * gammas(:, j, t)'), 2)', N, 1);
    d(:, :, j) = d(:, :, j) ./ repmat(sum(d(:, :, j), 2), 1, N);
%     if rank(d(:, :, j)) < N
%         [j, rank(d(:, :, j))]
%     end % if
end

% R_jn = bsxfun(@rdivide, w_nj .* L_nj, beta')';
% [mean(R_jn(:)), range(R_jn(:))]

% 
% 
% dd = reshape(permute(d, [3, 2, 1]), [N * J, N]);
% 
% for ii = 1:100
% E_jn = (gammas(:, :, t) + alpha(:, t) * beta') * R_jn  - alpha(:, t) * S(:, t)'; % * P_n_new(i_base);
% R_jn_new = reshape(sum(repmat(E_jn, [N, 1]) .* dd, 2), [J, N]);
% R_jn = R_jn_new;
% % [mean(R_jn(:)), range(R_jn(:))]
% end %for ii
% [mean(R_jn(:)), range(R_jn(:))]


B_d = [];

for n = 1:N
	temp = [];
	for m = 1:N
		d_vec = sparse(squeeze(d(m, n, :)));
		temp = [temp, diag(d_vec)];
	end % for m
	B_d = [B_d; temp];
end % for n


A = B_d * (B_gamma + D_alpha * B_beta) - eye(N * J);
b = B_d * D_alpha * S_full;
% b = zeros(size(S_full));

pinv_A = real(pinv(A));
x0 = pinv_A * b;
% e = A * x0 - b;

W = eye(N * J) - pinv_A * A;
% w1 = W * ones([(N * J) 1]);

% i0 = floor(-min(x0(x0 < 0) ./ sum(W(x0 < 0, :), 2))) + 1;

% aa = []; bb = [];
% for ii = 0:1000
% R = x0 + ii * W * w1;
% 
% R2 = x0 + 362.7090 * ii * W(:, 1);
% aa = [aa, sum(R)];
% bb = [bb, sum(R2)];
% end
% 
% plot([aa; bb]')

% R = x0 + i0 * w1;

ss = va_to_fit - sum(x0 .* beta_full);

scale_factor = ss / sum(W(:, 1) .* beta_full);

R = x0 + scale_factor * W(:, 1);

R_jn_new = reshape(R, [J, N]);

% va_nj_new = bsxfun(@times, R_jn_new', beta'); 

w_nj_new = real(bsxfun(@times, R_jn_new' ./ L_nj, beta'));


end
