function [w_nj_new, P_nj_new, price_iterations, P_n_new] = wage_update(w_nj, L_nj, z_nj, P_nj, t, va_to_fit, p_to_fit, S_t)
% This function is one iteration of the wage loop.                    (w_nj, L_nj, z_nj, P_nj, t, va_to_fit, p_to_fit, S_t)
% First it calculates prices corresponding to current wages,
% then calculates new wages.

global c rho nu betas theta kappa gammas 

assert_all(w_nj>0);
assert_all(L_nj>0);
assert_all(z_nj>0);
assert_all(P_nj>0);
assert_all(va_to_fit>0);
assert_all(p_to_fit>0);
assert_all(betas>0);
assert_all(kappa>=0);

i_base = c.i_base;
[N, J] = size(L_nj);

assert_all(size(S_t)==[N, 1]);
assert_all(size(gammas(:,:,t))==[J, J]);
assert_all(size(betas(:,t))==[J, 1]);
assert_all(abs(squeeze(sum(gammas, 1)) + betas - 1) < c.numerical_zero);

beta_full = repmat(betas(:,t), [N, 1]);

% Compute D to be used in the price equations. See the definition of D in the notes.
% D depends only on current sector specific wages and aggregate labor
D = bsxfun(@times, permute((bsxfun(@power, bsxfun(@times, sum(L_nj, 2), w_nj)', -betas(:,t)) .^ theta) .* z_nj', [3 2 1]), kappa(:, :, :, t).^theta);

assert_all(D>0);

% Get sectoral prices. Note that this step depends on current wages through D.
[P_nj_new, price_iterations] = get_prices(P_nj, D, t);

% check
% [mean(P_nj_new(:)), range(P_nj_new(:))]

% Calculate aggregate prices
P_n_new = aggregate_prices(P_nj_new, rho, nu(:,:,t));

% check
% [mean(P_n_new), range(P_n_new)]

% preallocate and compute matrix for d values to be used in the wage equations
d = zeros(N, N, J); 
for j = 1:J
    d(:, :, j) = D(:, :, j) .* repmat(prod(bsxfun(@power, P_nj_new, - theta * gammas(:, j, t)'), 2)', N, 1);
    d(:, :, j) = d(:, :, j) ./ repmat(sum(d(:, :, j), 2), 1, N);
%     if rank(d(:, :, j)) < N
%         [j, rank(d(:, :, j))]
%     end % if
end

%% FIXME: expanding a matrix on the fly is the slowest possible
B_d = [];
for n = 1:N
	temp = [];
	for m = 1:N
		d_vec = sparse(squeeze(d(m, n, :)));
		temp = [temp, diag(d_vec)];
	end % for m
	B_d = [B_d; temp];
end % for n

expenditure_shares = get_expenditure_shares(P_nj_new, rho, nu(:,:,t));
assert_all(size(expenditure_shares)==[N, J]);

%% FIXME: D_alpha depends on new prices
magic_A = zeros(N*J, N*J);
surplus_demand = zeros(N*J, 1);
for n = 1:N
	start_index = (n-1)*J+1;
	end_index = n*J;

	magic_A(start_index:end_index, start_index:end_index) = gammas(:,:,t) + expenditure_shares(n,:)'*betas(:,t)';
	surplus_demand(start_index:end_index) = S_t(n) * expenditure_shares(n,:)';
end


A = B_d * magic_A - eye(N * J);
b = B_d * surplus_demand;

options = optimoptions('lsqlin','Algorithm','active-set', 'Display', 'off');
R = lsqlin(A, b, [], [], beta_full', va_to_fit, [], [], [], options);

R_jn_new = reshape(R, [J, N]);

w_nj_new = real(bsxfun(@times, R_jn_new' ./ L_nj, betas(:,t)'));


end
