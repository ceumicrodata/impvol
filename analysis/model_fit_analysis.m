clear all

load('results/data_rgdp_and_volatility.mat')
load('results/equilibrium_baseline_theta_4.mat')

P_nt = equilibrium.P_nt;
w_nt = equilibrium.w_nt;
L_nt = equilibrium.L_nt;
P_njt = equilibrium.P_njt;
w_njt = equilibrium.w_njt;
L_njt = equilibrium.L_njt;
va_njt = w_njt .* L_njt;
va_nt = w_nt .* L_nt;

d_P_nt = deflator';
d_va_nt = va_total';
d_va_njt = va;
d_P_us_jt = p_sectoral_base';
d_P_njt = p_sectoral;

rgdp_nt = va_nt ./ P_nt;
d_rgdp_nt = d_va_nt ./ d_P_nt;

[N, J, T] = size(P_njt);



%% real GDP all countries
ROWS = 5;
COLS = 5;
for n = 1:N
    fig_number = floor((n - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(d_rgdp_nt(n, :))';
    series2 = log(rgdp_nt(n, :))';
    plot([series1, series2], 'LineWidth', 2)
    title([country_names{n}, ' Log Real GDP'])
    xlim([1, 36])
    ylim([7, 15])
    grid on
    set(gca, 'GridLineStyle', '-');
    if n == 1
        legend('Data', 'Model', 'Location', 'NW')
    end
    set(gca,'XTick', 1:8:36)
    set(gca,'XTickLabel', 1972:8:2007)
end
set(pic_id, 'Units', 'normalized', 'Position', [0, 0, 1, 1])
set(pic_id, 'PaperPositionMode','auto')
print(pic_id, 'realgdp', '-dpng', '-r0')
close(pic_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% deflator all countries
ROWS = 5;
COLS = 5;
for n = 1:N
    fig_number = floor((n - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(d_P_nt(n, :))';
    series2 = log(P_nt(n, :))';
    plot([series1, series2], 'LineWidth', 2)
    title([country_names{n}, ' Log Real GDP'])
    xlim([1, 36])
    ylim([1, 10])
    grid on
    set(gca, 'GridLineStyle', '-');
    if n == 1
        legend('Data', 'Model', 'Location', 'NW')
    end
    set(gca,'XTick', 1:8:36)
    set(gca,'XTickLabel', 1972:8:2007)
end
set(pic_id, 'Units', 'normalized', 'Position', [0, 0, 1, 1])
set(pic_id, 'PaperPositionMode','auto')
print(pic_id, 'deflator', '-dpng', '-r0')
close(pic_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% sectoral prices US
ROWS = 4;
COLS = 6;
for j = 1:J
    fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(j - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(d_P_us_jt(j, :))';
    series2 = log(squeeze(P_njt(25, j, :)));
    plot([series1, series2], 'LineWidth', 2)
    title(['Log Price, sec ', num2str(j)])
    xlim([1, 36])
    ylim([-1, 8])
    grid on
    set(gca, 'GridLineStyle', '-');
    if j == 1
        legend('Data', 'Model', 'Location', 'NW')
    end
    set(gca,'XTick', 1:8:36)
    set(gca,'XTickLabel', 1972:8:2007)
end
set(pic_id, 'Units', 'normalized', 'Position', [0, 0, 1, 1])
set(pic_id, 'PaperPositionMode','auto')
print(pic_id, 'usprices', '-dpng', '-r0')
close(pic_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% sectoral prices Mexico
ROWS = 4;
COLS = 6;
for j = 1:J
    fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(j - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(squeeze(d_P_njt(16, j, :)));
    series2 = log(squeeze(P_njt(16, j, :)));
    plot([series1, series2], 'LineWidth', 2)
    title(['Log Price, sec ', num2str(j)])
    xlim([1, 36])
    ylim([-1, 8])
    grid on
    set(gca, 'GridLineStyle', '-');
    if j == 1
        legend('Data', 'Model', 'Location', 'NW')
    end
    set(gca,'XTick', 1:8:36)
    set(gca,'XTickLabel', 1972:8:2007)
end
set(pic_id, 'Units', 'normalized', 'Position', [0, 0, 1, 1])
set(pic_id, 'PaperPositionMode','auto')
print(pic_id, 'mexicoprices', '-dpng', '-r0')
close(pic_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %% sectoral value added US
% ROWS = 4;
% COLS = 6;
% for j = 1:J
%     fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
%     pic_id = figure(fig_number);
%     hold on
%     thisplot = mod(j - 1, ROWS * COLS) + 1;
%     subplot(ROWS, COLS, thisplot)
%     series1 = log(squeeze(d_va_njt(25, j, :)));
%     series2 = log(squeeze(va_njt(25, j, :)));
%     plot([series1, series2], 'LineWidth', 2)
%     title(['Log Value Added, sec ', num2str(j)])
%     xlim([1, 36])
%     ylim([7, 21])
%     grid on
%     set(gca, 'GridLineStyle', '-');
%     if j == 1
%         legend('Data', 'Model', 'Location', 'NW')
%     end
%     set(gca,'XTick', 1:8:36)
%     set(gca,'XTickLabel', 1972:8:2007)
% end
% set(pic_id, 'Units', 'normalized', 'Position', [0, 0, 1, 1])
% set(pic_id, 'PaperPositionMode','auto')
% print(pic_id, 'usva', '-dpng', '-r0')
% close(pic_id)
