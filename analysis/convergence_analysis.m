load('iterations.mat')

[N, J, T, iterations] = size(L_njt_iterations);
iterations = iterations - 1;

L_nt_full = permute(repmat(L_nt, [1, 1, J, iterations + 1]), [1 3 2 4]);
L_njt_shares_iterations = L_njt_iterations ./ L_nt_full;

diffs = abs(diff(L_njt_iterations, 1, 4));
reldiffs = abs(diff(L_njt_shares_iterations, 1, 4));

nn = 25;
jj = 24;

figure(1)
plot(squeeze(L_njt_shares_iterations(nn, jj, :, :)))
figure(2)
plot(squeeze(reldiffs(nn, jj, 1, :)))

