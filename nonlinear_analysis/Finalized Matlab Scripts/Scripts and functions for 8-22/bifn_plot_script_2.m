% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 3; 

grid_row = 3;
grid_col = 5;

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

a_val_1 = 0;
a_val_2 = 1/12;
a_val_3 = 1/6;

% fully recessive case ----------------------------------------------------

k_val = 1; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 1;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 6;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 11;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% partially recessive case ------------------------------------------------

k_val = .75; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 2;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 7;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 12;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% additive case -----------------------------------------------------------

k_val = .5; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 3;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 8;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 13;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% partially dominant case -----------------------------------------------------

k_val = .25; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 4;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 9;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 14;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% fully dominant case -----------------------------------------------------

k_val = 0; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 5;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 10;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 15;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

savefig("15_bifn_diagram_v1.fig")



% Functions: 


%plotting functions -------------------------------------------------------

function ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if h_val <= 2/3

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

        if h_val == 0
            legend('0 HEs', 'α=0', 'Diploids')
        end

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



function autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [autos1_q1, autos1_s1, autos1_q2, autos1_s2, autos1_q3, autos1_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_1);

    [autos2_q1, autos2_s1, autos2_q2, autos2_s2, autos2_q3, autos2_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_2);
  
    [autos3_q1, autos3_s1, autos3_q2, autos3_s2, autos3_q3, autos3_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_3);

    if h_val <= 2/3

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


