function plot_trade_costs(kappa)

[N, ~, J, T] = size(kappa);

kappa_agr_1972 = squeeze(kappa(:, :, 1, 1));
temp = triu(kappa_agr_1972, 1) - tril(ones(size(kappa_agr_1972)));
kappa_agr_1972_v = temp(temp >= 0);

kappa_agr_2007 = squeeze(kappa(:, :, 1, T));
temp = triu(kappa_agr_2007, 1) - tril(ones(size(kappa_agr_2007)));
kappa_agr_2007_v = temp(temp >= 0);

figure1 = figure;
subplot(1, 2, 1)
hist(kappa_agr_1972_v, 0.002:0.004:0.3)
xlim([0 0.2])
ylim([0 60])
title('1972')
subplot(1, 2, 2)
hist(kappa_agr_2007_v, 0.002:0.004:0.3)
xlim([0 0.2])
ylim([0 60])
title('2007')
set(figure1, 'Units','normalized','position',[.2 .2 .4 .3])
set(figure1,'PaperPositionMode','auto')
print(figure1, '-dpng', '-r600', 'results\\fig4.png') 
close(figure1)

kappa_manuf_1972_v = zeros(N * (N - 1) / 2, J - 2);
kappa_manuf_2007_v = zeros(N * (N - 1) / 2, J - 2);

for j = 2:(J - 1)
    kappa_manuf_j_1972 = squeeze(kappa(:, :, j, 1));
    temp = triu(kappa_manuf_j_1972, 1) - tril(ones(size(kappa_manuf_j_1972)));
    kappa_manuf_1972_v(:, j - 1) = temp(temp >= 0);
    
    kappa_manuf_j_2007 = squeeze(kappa(:, :, j, T));
    temp = triu(kappa_manuf_j_2007, 1) - tril(ones(size(kappa_manuf_j_2007)));
    kappa_manuf_2007_v(:, j - 1) = temp(temp >= 0);
end % j

figure2 = figure;
subplot(1, 2, 1)
hist(kappa_manuf_1972_v(:), 0.01:0.02:1)
xlim([0 1])
ylim([0 2600])
title('1972')
subplot(1, 2, 2)
hist(kappa_manuf_2007_v(:), 0.01:0.02:1)
xlim([0 1])
ylim([0 2600])
title('2007')
set(figure2, 'Units','normalized','position',[.2 .2 .4 .3])
set(figure2,'PaperPositionMode','auto')
print(figure2, '-dpng', '-r600', 'results\\fig3.png') 
close(figure2)

end


