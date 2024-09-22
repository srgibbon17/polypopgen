function autos_F_stat_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [autos1_F1, autos1_s1, autos1_F2, autos1_s2, autos1_F3, autos1_s3] = auto_F_stat_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_1);

    [autos2_F1, autos2_s1, autos2_F2, autos2_s2, autos2_F3, autos2_s3] = auto_F_stat_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_2);
  
    [autos3_F1, autos3_s1, autos3_F2, autos3_s2, autos3_F3, autos3_s3] = auto_F_stat_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_3);

    if isempty(autos1_F3) == 1 && isempty(autos2_F3) == 1 && isempty(autos3_F3) == 1

        autos_1_q_values = cat(2, autos1_F1, autos1_F2);
        autos_1_s_values = cat(2, autos1_s1, autos1_s2);

        autos_2_q_values = cat(2, autos2_F1, autos2_F2);
        autos_2_s_values = cat(2, autos2_s1, autos2_s2);

        autos_3_q_values = cat(2, autos3_F1, autos3_F2);
        autos_3_s_values = cat(2, autos3_s1, autos3_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(autos_1_s_values, autos_1_q_values, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.25)
        hold on
        plot(autos_2_s_values, autos_2_q_values, 'Color', 'k', 'LineWidth', 1.25)
        plot(autos_3_s_values, autos_3_q_values, 'Color', "#0072BD", 'LineWidth', 1.25)

        if h_val == 0
            legend({'α=0', 'α=1/12', 'α=1/6'}, 'FontSize', 16)
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(autos1_s1, autos1_F1, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.25)
        hold on
        plot(autos1_s2, autos1_F2, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.25)
        plot(autos1_s3, autos1_F3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--', 'LineWidth', 1.25)

        plot(autos2_s1, autos2_F1, 'Color', 'k', 'LineWidth', 1.25)
        plot(autos2_s2, autos2_F2, 'Color', 'k', 'LineWidth', 1.25)
        plot(autos2_s3, autos2_F3, 'Color', 'k', 'LineStyle','--', 'LineWidth', 1.25)

        plot(autos3_s1, autos3_F1, 'Color', "#0072BD", 'LineWidth', 1.25)
        plot(autos3_s2, autos3_F2, 'Color', "#0072BD", 'LineWidth', 1.25)
        plot(autos3_s3, autos3_F3, 'Color', "#0072BD", 'LineStyle','--', 'LineWidth', 1.25)

    end
    
    xscale log

    if grid_pos == 1
        ylabel('F-statistic', 'FontSize', 24)
    end

    if grid_pos == 3
        xlabel('s (selection coefficient)', 'FontSize', 24)
    end

    ylim([-1, 1])

    disp(grid_pos)

end

