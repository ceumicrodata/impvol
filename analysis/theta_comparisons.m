clear all

load('results\data_rgdp_and_volatility.mat', 'realGdpTotal')
load('algorithm_input\data_theta_2.mat', 'baseline')
countryNames = baseline.countryNames;

rGDP_data = realGdpTotal;
[N, T] = size(rGDP_data);

rGDP = zeros(N, T, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          THETA 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results\equilibrium_baseline_theta_2.mat')
w_nt = equilibrium(1).w_nt;
L_nt = equilibrium(1).L_nt;
P_nt = equilibrium(1).P_nt;
total_value_added = w_nt .* L_nt;
rGDP(:, :, 1) = total_value_added ./ P_nt;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          THETA 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results\equilibrium_baseline_theta_4.mat')
w_nt = equilibrium(1).w_nt;
L_nt = equilibrium(1).L_nt;
P_nt = equilibrium(1).P_nt;
total_value_added = w_nt .* L_nt;
rGDP(:, :, 2) = total_value_added ./ P_nt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          THETA 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('results\equilibrium_baseline_theta_8.mat')
w_nt = equilibrium(1).w_nt;
L_nt = equilibrium(1).L_nt;
P_nt = equilibrium(1).P_nt;
total_value_added = w_nt .* L_nt;
rGDP(:, :, 3) = total_value_added ./ P_nt;



%%%%%%%%%%%%%    PLOT    %%%%%%%%%%%%%%%%%%%%%%

ROWS = 5;
COLS = 5;
figure()
for n = 1:N
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(rGDP_data(n, :))';
    series2 = squeeze(log(rGDP(n, :, :)));
    plot([series1, series2], 'LineWidth', 1.5)
    title([countryNames{n}, ' Log rGDP'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    if n == 1
        legend('Data', 'Theta 2', 'Theta 4', 'Theta 8', 'Location', 'SouthEast')
    end
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
end