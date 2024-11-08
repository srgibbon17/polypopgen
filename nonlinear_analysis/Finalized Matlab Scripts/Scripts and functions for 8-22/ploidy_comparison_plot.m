function ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if isempty(diploids_q3) == 1 && isempty(autos_q3) == 1 && isempty(HE0_q3) == 1 

        HE0_q_values = cat(2, HE0_q1, HE0_q2);
        HE0_s_values = cat(2, HE0_s1, HE0_s2);

        autos_q_values = cat(2, autos_q1, autos_q2);
        autos_s_values = cat(2, autos_s1, autos_s2);

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s_values, HE0_q_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(autos_s_values, autos_q_values, 'Color', [0.8500 0.3250 0.0980])
        plot(diploids_s_values, diploids_q_values, 'Color', 'k')

       
        legend('0 HEs', 'α=0', 'Diploids')
        

    elseif isempty(diploids_q3) == 1 && isempty(autos_q3) == 0 && isempty(HE0_q3) == 0 

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(diploids_s_values, diploids_q_values, 'Color', 'k')
        hold on
        plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(autos_s1, autos_q1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s2, autos_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s3, autos_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(autos_s1, autos_q1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s2, autos_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s3, autos_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(diploids_s1, diploids_q1, 'Color','k')
        plot(diploids_s2, diploids_q2, 'Color','k')
        plot(diploids_s3, diploids_q3, 'Color','k', 'LineStyle','--')

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