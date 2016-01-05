load('D:\peter\impvol\results_backup\new\equilibrium_theta_4.mat')
load('D:\peter\impvol\results_backup\new\equilibrium_6_theta_4.mat', 'counterfactual_equilibrium')
d_benchmark = equilibrium(1, 1).d;
d_free_trade = equilibrium(1, 7).d;
load('D:\peter\impvol\results_backup\equilibrium_3_theta_4.mat')
d_1927 = counterfactual_equilibrium.d;


dd = d_data;
for j=1:23
    for t=1:36
        D(:, j, t) = diag(dd(:, :, j, t));
    end
end
disp('DATA')
mean(D(:))

dd = d_benchmark;
for j=1:23
    for t=1:36
        D(:, j, t) = diag(dd(:, :, j, t));
    end
end
disp('BENCHMARK')
mean(D(:))

dd = d_free_trade;
for j=1:23
    for t=1:36
        D(:, j, t) = diag(dd(:, :, j, t));
    end
end
disp('FREE TRADE')
mean(D(:))

dd = d_1972;
for j=1:23
    for t=1:36
        D(:, j, t) = diag(dd(:, :, j, t));
    end
end
disp('1972')
mean(D(:))


hist(
