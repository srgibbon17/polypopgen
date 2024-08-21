% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 100; 

grid_row = 1;
grid_col = 3;

% recessive case:

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 1; % double reduction rate

h_val = 0; % diploid dominance coefficient

h1_val = 0; % simplex dominance coefficient
h2_val = 0; % duplex dominance coefficient
h3_val = 0; % triplex dominance coefficient

grid_pos = 1;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)


% additive case:

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 1; % double reduction rate

h_val = .5; % diploid dominance coefficient

h1_val = .25; % simplex dominance coefficient
h2_val = .5; % duplex dominance coefficient
h3_val = .75; % triplex dominance coefficient

grid_pos = 2;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% additive case:

s_val_range = logspace(-8, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 1; % double reduction rate

h_val = 1; % diploid dominance coefficient

h1_val = 1; % simplex dominance coefficient
h2_val = 1; % duplex dominance coefficient
h3_val = 1; % triplex dominance coefficient

grid_pos = 3;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)


function ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    %[HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE2_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if h_val <= 2/3

        %HE0_q_values = cat(2, HE0_q1, HE0_q2);
        %HE0_s_values = cat(2, HE0_s1, HE0_s2);

        autos_q_values = cat(2, autos_q1, autos_q2);
        autos_s_values = cat(2, autos_s1, autos_s2);

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        %plot(HE0_s_values, HE0_q_values, 'Color', [0 0.4470 0.7410])
        
        plot(autos_s_values, autos_q_values, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(diploids_s_values, diploids_q_values, 'Color', 'k')

        if h_val == 0
            legend('α=0', 'Diploids')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        %plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        
        %plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        %plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(autos_s1, autos_q1, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(autos_s2, autos_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s3, autos_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(diploids_s1, diploids_q1, 'Color','k')
        plot(diploids_s2, diploids_q2, 'Color','k')
        plot(diploids_s3, diploids_q3, 'Color','k', 'LineStyle','--')

    end


    if h_val == 0
        title('h=0 (recessive)')
    elseif h_val == .5
        title('h=.5 (additive)')
    elseif h_val == 1
        title('h=1 (dominant)')
    end
    
    xscale log

    if grid_pos == 1 
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 2
        xlabel('s (selection coefficient)')
    end
    
    ylim([0, 1])

end
