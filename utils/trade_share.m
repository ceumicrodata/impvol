
clear all

%% Baseline
m0 = 'table_th_4_lac_inf_baseline';
e0 = get_eq_trade(m0);

%% Kappa 1972 counterfactual 
m1 = 'table_th_4_lac_inf_kappa1972';
e1 = get_eq_trade(m1);

%% check (1-beta)/beta
gamma = e0.gamma(:, :, 1);
beta = e0.beta;

%% Actual trade volumes
e0.share_nt = e0.I_nt ./ e0.va_nt;
e1.share_nt = e1.I_nt ./ e1.va_nt;

figure()
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.share_nt(n, :); e1.share_nt(n, :)]')
    title(e0.names{n})
    legend('Baseline', '1972', 'Location', 'NW')
end