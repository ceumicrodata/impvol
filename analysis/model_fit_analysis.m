clear all
% load('results\data)
load('results\equilibrium_baseline_theta_4.mat')
% load('data\sandbox_algorithm_input\data_theta_2.mat')

% load('results\data_rgdp_and_volatility.mat')
% load('results\equilibrium_baseline_theta_2.mat')
load('algorithm_input\data_theta_4.mat')

country_names = baseline.country_names;

w_njt = equilibrium(1).w_njt;
L_njt = equilibrium(1).L_njt;

w_nt = equilibrium(1).w_nt;
L_nt = equilibrium(1).L_nt;

P_nt = equilibrium(1).P_nt;
P_njt = equilibrium(1).P_njt;

[N, J, T] = size(L_njt);
%i_base = N;

sectoral_value_added = w_njt .* L_njt;
total_value_added = w_nt .* L_nt;

% aa = squeeze(sum(sectoral_value_added, 2));
% max(abs(aa(:) - total_value_added(:)));

rGDP = total_value_added ./ P_nt;
rGDP_sectoral = sectoral_value_added ./ permute(repmat(P_nt, [1 1 J]), [1 3 2]);

rGDP_data = real_gdp_total;
rGDP_sectoral_data = real_dgp_sectoral; 


% Real GDP
ROWS = 5;
COLS = 5;
for n = 1:N
    fig_number = floor((n - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series2 = log(rGDP(n, :))';
    series1 = log(rGDP_data(n, :))';
%     plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, series2], 'LineWidth', 2.5)
%     plot([series2], 'LineWidth', 2.5)
    title([country_names{n}, ' Log Real GDP'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    %legend('Data', 'Model', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0) || (n == 25)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
            export_fig('results\\all_countries_log_rGDP_fe_alpha.pdf', '-append')
        else
            export_fig('results\\all_countries_log_rGDP_fe_alpha.pdf')
        end % fig_number
        close(pic_id)
    end % if
end

% Value Added
ROWS = 2;
COLS = 3;
for n = 1:N
    fig_number = floor((n - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(total_value_added(n, :)');
    series2 = log(va_total(:, n));
%     plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, series2], 'LineWidth', 2.5)
    %plot(zscore([series1; series2]'), 'LineWidth', 2.5)
    title([country_names{n}, ' Log Value Added, baseline model'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    legend('Model', 'Data', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0) || (n == 25)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
            %export_fig('sandbox_results\\all_countries_VA_baseline.pdf', '-append')
        else
            %export_fig('sandbox_results\\all_countries_VA_baseline.pdf')
        end % fig_number
        %close(pic_id)
    end % if
end


% GDP Deflator
ROWS = 5;
COLS = 5;
for n = 1:N
    fig_number = floor((n - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(n - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(P_nt(n, :)');
    series2 = log(deflator(:, n));
%     plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, series2], 'LineWidth', 2.5)
    %plot(zscore([series1; series2]'), 'LineWidth', 2.5)
    title([country_names{n}, ' Log GDP Deflator'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    legend('Model', 'Data', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0) || (n == 25)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
%             export_fig('results\\all_countries_log_GDP_deflator_fe_alpha.pdf', '-append')
        else
%             export_fig('results\\all_countries_log_GDP_deflator_fe_alpha.pdf')
        end % fig_number
%         close(pic_id)
    end % if
end

% Sectoral Value Added in a given Country
ROWS = 4;
COLS = 6;
for n = 1:N
    figure(n)
    for j = 1:J
%         fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
        
%         pic_id = figure(fig_number);
        hold on
        thisplot = mod(j - 1, ROWS * COLS) + 1;
        subplot(ROWS, COLS, thisplot)
        series1 = log(squeeze(sectoral_value_added(n, j, :)));
        series2 = log(squeeze(va(n, j, :)));
        %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
        plot([series1, series2], 'LineWidth', 2.5)
        title([country_names{n}, ', Sec ', num2str(j)])
        xlim([1, 36])
        grid on
        set(gca, 'GridLineStyle', '-');
        legend('Model', 'Data', 'Location', 'SouthEast')
        set(gca,'XTick', 1:5:36)
        set(gca,'XTickLabel', 1972:5:2007)
        if (mod(thisplot, ROWS * COLS) == 0)
            set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
            if fig_number > 1
%                 export_fig(['results\\', country_names{n}, '_sectoral_log_va_fe_alpha.pdf'], '-append')
            else
%                 export_fig(['results\\', country_names{n}, '_sectoral_log_va_fe_alpha.pdf'])
            end % fig_number
%             close(pic_id)
        end % if
    end % j
    hold off
end % n

% Sectoral Real GDP in a given Country
ROWS = 4;
COLS = 6;
for n = 1:N
    for j = 1:J
        fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
        pic_id = figure(fig_number);
        hold on
        thisplot = mod(j - 1, ROWS * COLS) + 1;
        subplot(ROWS, COLS, thisplot)
        series1 = log(squeeze(rGDP_sectoral(n, j, :)));
        series2 = log(squeeze(rGDP_sectoral_data(n, j, :)));
        %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
        plot([series1, series2], 'LineWidth', 2.5)
        title([country_names{n}, ', Sec ', num2str(j)])
        xlim([1, 36])
        grid on
        set(gca, 'GridLineStyle', '-');
        legend('Model', 'Data', 'Location', 'SouthEast')
        set(gca,'XTick', 1:5:36)
        set(gca,'XTickLabel', 1972:5:2007)
        if (mod(thisplot, ROWS * COLS) == 0)
            set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
            if fig_number > 1
                export_fig(['results\\', country_names{n}, '_sectoral_rGDP_fe_alpha.pdf'], '-append')
            else
                export_fig(['results\\', country_names{n}, '_sectoral_rGDP_fe_alpha.pdf'])
            end % fig_number
            close(pic_id)
        end % if
    end % j
end % n

% Sectoral Prices in the US
ROWS = 2;
COLS = 3;
for j = 1:J
    fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(j - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(squeeze(P_njt(iUs, j, :)));
    series2 = log(p_sectoralUs(:, j));   
    %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, series2], 'LineWidth', 2.5)
    title(['UK Prices in Sector ', num2str(j), ', baseline model'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    legend('Model', 'Data', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
            export_fig('sandbox_results\\UK_sectoral_prices_baseline.pdf', '-append')
        else
            export_fig('sandbox_results\\UK_sectoral_prices_baseline.pdf')
        end % fig_number
        close(pic_id)
    end % if
end

% Sectoral Value Added in the US
ROWS = 2;
COLS = 3;
for j = 1:J
    fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(j - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = log(squeeze(sectoral_value_added(iUs, j, :)));
    series2 = log(squeeze(va(iUs, j, :)));  
    %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, series2], 'LineWidth', 2.5)
    title(['US Value Added in Sector ', num2str(j), ', baseline model'])
    xlim([1, 36])
    grid on
    set(gca, 'GridLineStyle', '-');
    legend('Model', 'Data', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
            export_fig('sandbox_results\\US_sectoral_VA_baseline.pdf', '-append')
        else
            export_fig('sandbox_results\\US_sectoral_VA_baseline.pdf')
        end % fig_number
        close(pic_id)
    end % if
end




% Sectoral Productivity in a given Country
ROWS = 2;
COLS = 3;
for n = 1:N
    for j = 1:J
        fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
        pic_id = figure(fig_number);
        hold on
        thisplot = mod(j - 1, ROWS * COLS) + 1;
        subplot(ROWS, COLS, thisplot)
        series1 = log(squeeze(z_us(n, j, :)));
        series2 = log(squeeze(z_dk(n, j, :)));
        %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
        plot([series1, series2], 'LineWidth', 2.5)
        title([country_names{n}, ', Sector ', num2str(j), ', Log Productivity'])
        xlim([1, 36])
        grid on
        set(gca, 'GridLineStyle', '-');
        legend('US', 'Mean d*kappa', 'Location', 'SouthEast')
        set(gca,'XTick', 1:5:36)
        set(gca,'XTickLabel', 1972:5:2007)
        if (mod(thisplot, ROWS * COLS) == 0)
            set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
            if fig_number > 1
                export_fig(['sandbox_results\\', country_names{n}, '_sectoral_productivity.pdf'], '-append')
            else
                export_fig(['sandbox_results\\', country_names{n}, '_sectoral_productivity.pdf'])
            end % fig_number
            close(pic_id)
        end % if
    end % j
end % n


%%%%%%%     Sectors     %%%%%%%%%%%%

va_share_model = ...
    bsxfun(@rdivide, sectoral_value_added, sum(sectoral_value_added, 2));
va_share_data = bsxfun(@rdivide, va, sum(va, 2));

% Sectoral Prices in the US
ROWS = 4;
COLS = 6;
for j = 1:J
    fig_number = floor((j - 1) / (ROWS * COLS)) + 1;
    pic_id = figure(fig_number);
    hold on
    thisplot = mod(j - 1, ROWS * COLS) + 1;
    subplot(ROWS, COLS, thisplot)
    series1 = squeeze(...
        mean(log(va_share_model(:, j, :) ./ va_share_data(:, j, :)), 1));
%     series2 = log(p_sectoralUs(:, j));   
    %plot([series1 - mean(series1), series2 - mean(series2)], 'LineWidth', 2.5)
    plot([series1, zeros(36, 1)], 'LineWidth', 2.5)
    title(['Diff in VA share, Sector ', num2str(j),])
    xlim([1, 36])
    ylim([-0.5, 1.75])
    grid on
    set(gca, 'GridLineStyle', '-');
%     legend('Model', 'Data', 'Location', 'SouthEast')
    set(gca,'XTick', 1:5:36)
    set(gca,'XTickLabel', 1972:5:2007)
    if (mod(thisplot, ROWS * COLS) == 0)
        set(pic_id, 'Units','normalized','position',[.1 .1 .9 .9])
        if fig_number > 1
            export_fig('results\\sectoral_va_share.pdf', '-append')
        else
            export_fig('results\\sectoral_va_share.pdf')
        end % fig_number
        close(pic_id)
    end % if
end
