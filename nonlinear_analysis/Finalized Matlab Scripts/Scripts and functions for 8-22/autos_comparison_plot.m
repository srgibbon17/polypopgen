function autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [autos1_q1, autos1_s1, autos1_q2, autos1_s2, autos1_q3, autos1_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_1);

    [autos2_q1, autos2_s1, autos2_q2, autos2_s2, autos2_q3, autos2_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_2);
  
    [autos3_q1, autos3_s1, autos3_q2, autos3_s2, autos3_q3, autos3_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_3);

    if isempty(autos1_q3) == 1 && isempty(autos2_q3) == 1 && isempty(autos3_q3) == 1 

        autos_1_q_values = cat(2, autos1_q1, autos1_q2);
        autos_1_s_values = cat(2, autos1_s1, autos1_s2);

        autos_2_q_values = cat(2, autos2_q1, autos2_q2);
        autos_2_s_values = cat(2, autos2_s1, autos2_s2);

        autos_3_q_values = cat(2, autos3_q1, autos3_q2);
        autos_3_s_values = cat(2, autos3_s1, autos3_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(autos_1_s_values, autos_1_q_values, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(autos_2_s_values, autos_2_q_values, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos_3_s_values, autos_3_q_values, 'Color', [0.8500 0.3250 0.0980 .4])

        if h_val == 0
            legend('α=0', 'α=1/12', 'α=1/6')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(autos1_s1, autos1_q1, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(autos1_s2, autos1_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s3, autos1_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(autos2_s1, autos2_q1, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos2_s2, autos2_q2, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos2_s3, autos2_q3, 'Color', [0.8500 0.3250 0.0980 .7], 'LineStyle','--')

        plot(autos3_s1, autos3_q1, 'Color', [0.8500 0.3250 0.0980 .4])
        plot(autos3_s2, autos3_q2, 'Color', [0.8500 0.3250 0.0980 .4])
        plot(autos3_s3, autos3_q3, 'Color', [0.8500 0.3250 0.0980 .4], 'LineStyle','--')

    end
    
    xscale log

    if grid_pos == 6
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 13
        xlabel('s (selection coefficient)')
    end

    ylim([0, 1])

    disp(grid_pos)

end

