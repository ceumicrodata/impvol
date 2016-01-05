function [w_njt_new, w_nt_new, P_nt_new, P_njt_new] = fit_to_data(w_njt, w_nt, P_nt, P_njt)

global c

inFolder = c.resultsFolder;
iBase = c.iBase;

load([inFolder, 'data_rgdp_and_volatility.mat'], 'pBase')

[nCountries, nSectors, ~] = size(P_njt);


correction = repmat(pBase' ./ P_nt(iBase, :), [nCountries 1]);
correction_full = permute(repmat(correction, [1 1 nSectors]), [1 3 2]);

P_nt_new = P_nt .* correction;
w_nt_new = w_nt .* correction;

P_njt_new = P_njt .* correction_full;
w_njt_new = w_njt .* correction_full;

end