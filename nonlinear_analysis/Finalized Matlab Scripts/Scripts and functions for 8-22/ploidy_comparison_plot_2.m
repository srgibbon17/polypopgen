function ploidy_comparison_plot_2(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

    [allo_q1, allo_s1, allo_q2, allo_s2, allo_q3, allo_s3] = allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data_2(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if isempty(diploids_q3) == 1 && isempty(autos_q3) == 1 && isempty(allo_q3) == 1 

        allo_q_values = cat(2, allo_q1, allo_q2);
        allo_s_values = cat(2, allo_s1, allo_s2);

        autos_q_values = cat(2, autos_q1, autos_q2);
        autos_s_values = cat(2, autos_s1, autos_s2);

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(allo_s_values, allo_q_values, 'Color', allo_color, 'LineWidth', 1)
        hold on
        plot(autos_s_values, autos_q_values, 'Color', auto_color, 'LineWidth', 1)
        plot(diploids_s_values, diploids_q_values, 'Color', diploid_color, 'LineWidth', 1)

        if grid_pos == 6
            legend('Allos', 'Autos', 'Diploids')
        end

    elseif isempty(diploids_q3) == 1 && isempty(autos_q3) == 0 && isempty(allo_q3) == 0 

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(diploids_s_values, diploids_q_values, 'Color', diploid_color, 'LineWidth', 1)
        hold on
        plot(allo_s1, allo_q1, 'Color', allo_color, 'LineWidth', 1)
        plot(allo_s2, allo_q2, 'Color', allo_color, 'LineWidth', 1)
        plot(allo_s3, allo_q3, 'Color', allo_color, 'LineStyle','--', 'LineWidth', 1)

        plot(autos_s1, autos_q1, 'Color', auto_color, 'LineWidth', 1)
        plot(autos_s2, autos_q2, 'Color', auto_color, 'LineWidth', 1)
        plot(autos_s3, autos_q3, 'Color', auto_color, 'LineStyle','--', 'LineWidth', 1)

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(allo_s1, allo_q1, 'Color', allo_color, 'LineWidth', 1)
        hold on
        plot(allo_s2, allo_q2, 'Color', allo_color, 'LineWidth', 1)
        plot(allo_s3, allo_q3, 'Color', allo_color, 'LineStyle','--', 'LineWidth', 1)

        plot(autos_s1, autos_q1, 'Color', auto_color, 'LineWidth', 1)
        plot(autos_s2, autos_q2, 'Color', auto_color, 'LineWidth', 1)
        plot(autos_s3, autos_q3, 'Color', auto_color, 'LineStyle','--', 'LineWidth', 1)

        plot(diploids_s1, diploids_q1, 'Color', diploid_color, 'LineWidth', 1)
        plot(diploids_s2, diploids_q2, 'Color', diploid_color, 'LineWidth', 1)
        plot(diploids_s3, diploids_q3, 'Color', diploid_color, 'LineStyle','--', 'LineWidth', 1)

    end

    xscale log

    if grid_pos == 1 || grid_pos == 4 
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 5
        xlabel('s (selection coefficient)')
    end
    
    ylim([0, 1])

    disp(grid_pos)

end