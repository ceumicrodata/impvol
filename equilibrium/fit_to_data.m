function [w_njt_new, w_nt_new, P_nt_new, P_njt_new] = fit_to_data(w_njt, w_nt, P_nt, P_njt)

global c

in_folder = c.results_folder;
i_base = c.i_base;

load([in_folder, 'data_rgdp_and_volatility.mat'], 'p_base')

[n_countries, n_sectors, ~] = size(P_njt);


correction = repmat(p_base' ./ P_nt(i_base, :), [n_countries 1]);
correction_full = permute(repmat(correction, [1 1 n_sectors]), [1 3 2]);

P_nt_new = P_nt .* correction;
w_nt_new = w_nt .* correction;

P_njt_new = P_njt .* correction_full;
w_njt_new = w_njt .* correction_full;

end