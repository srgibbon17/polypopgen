function allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [HE1_q1, HE1_s1, HE1_q2, HE1_s2, HE1_q3, HE1_s3] = HE1_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [HE2_q1, HE2_s1, HE2_q2, HE2_s2, HE2_q3, HE2_s3] = HE2_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    if h_val <= 2/3

        HE0_q_values = cat(2, HE0_q1, HE0_q2);
        HE0_s_values = cat(2, HE0_s1, HE0_s2);

        HE1_q_values = cat(2, HE1_q1, HE1_q2);
        HE1_s_values = cat(2, HE1_s1, HE1_s2);

        HE2_q_values = cat(2, HE2_q1, HE2_q2);
        HE2_s_values = cat(2, HE2_s1, HE2_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s_values, HE0_q_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE1_s_values, HE1_q_values, 'Color', [0 0.4470 0.7410 .7])
        plot(HE2_s_values, HE2_q_values, 'Color', [0 0.4470 0.7410 .4])

        if h_val == 0
            legend('0 HEs', '1 HE', '2 HEs')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(HE1_s1, HE1_q1, 'Color', [0 0.4470 0.7410 .7])
        plot(HE1_s2, HE1_q2, 'Color', [0 0.4470 0.7410 .7])
        plot(HE1_s3, HE1_q3, 'Color', [0 0.4470 0.7410 .7], 'LineStyle','--')

        plot(HE2_s1, HE2_q1, 'Color', [0 0.4470 0.7410 .4])
        plot(HE2_s2, HE2_q2, 'Color', [0 0.4470 0.7410 .4])
        plot(HE2_s3, HE2_q3, 'Color', [0 0.4470 0.7410 .4], 'LineStyle','--')

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