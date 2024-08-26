% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 15; 

grid_row = 3;
grid_col = 3;

% constant parameters
mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

a_val_1 = 0;
a_val_2 = 1/4;

% recessive case:

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

h_val = 0; % diploid dominance coefficient

h1_val = 0; % simplex dominance coefficient
h2_val = 0; % duplex dominance coefficient
h3_val = 0; % triplex dominance coefficient

grid_pos = 1;

bifn_eigen_stiff_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% additive case:

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

h_val = .5; % diploid dominance coefficient

h1_val = .25; % simplex dominance coefficient
h2_val = .5; % duplex dominance coefficient
h3_val = .75; % triplex dominance coefficient

grid_pos = 2;

bifn_eigen_stiff_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% dominant case:

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

h_val = 1; % diploid dominance coefficient

h1_val = 1; % simplex dominance coefficient
h2_val = 1; % duplex dominance coefficient
h3_val = 1; % triplex dominance coefficient

grid_pos = 3;

bifn_eigen_stiff_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)


function bifn_eigen_stiff_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [allos1_q1, allos1_s1, allos1_eig1, allos1_stiff1, allos1_q2, allos1_s2, allos1_eig2, allos1_stiff2, allos1_q3, allos1_s3] = HE0_eigen_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [allos2_q1, allos2_s1, allos2_eig1, allos2_stiff1, allos2_q2, allos2_s2, allos2_eig2, allos2_stiff2, allos2_q3, allos2_s3] = HE2_eigen_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos1_q1, autos1_s1, autos1_eig1, autos1_stiff1, autos1_q2, autos1_s2, autos1_eig2, autos1_stiff2, autos1_q3, autos1_s3] = auto_eigen_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_1);

    [autos2_q1, autos2_s1, autos2_eig1, autos2_stiff1, autos2_q2, autos2_s2, autos2_eig2, autos2_stiff2, autos2_q3, autos2_s3] = auto_eigen_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_2);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if h_val <= 2/3

        allos1_q_values = cat(2, allos1_q1, allos1_q2);
        allos1_s_values = cat(2, allos1_s1, allos1_s2);
        allos1_eig_values = cat(2, allos1_eig1, allos1_eig2);
        allos1_stiff_values = cat(2, allos1_stiff1, allos1_stiff2);

        allos2_q_values = cat(2, allos2_q1, allos2_q2);
        allos2_s_values = cat(2, allos2_s1, allos2_s2);
        allos2_eig_values = cat(2, allos2_eig1, allos2_eig2);
        allos2_stiff_values = cat(2, allos2_stiff1, allos2_stiff2);

        autos1_q_values = cat(2, autos1_q1, autos1_q2);
        autos1_s_values = cat(2, autos1_s1, autos1_s2);
        autos1_eig_values = cat(2, autos1_eig1, autos1_eig2);
        autos1_stiff_values = cat(2, autos1_stiff1, autos1_stiff2);

        autos2_q_values = cat(2, autos2_q1, autos2_q2);
        autos2_s_values = cat(2, autos2_s1, autos2_s2);
        autos2_eig_values = cat(2, autos2_eig1, autos2_eig2);
        autos2_stiff_values = cat(2, autos2_stiff1, autos2_stiff2);

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        % bifn plot
        subplot(grid_row, grid_col, grid_pos)

        plot(allos1_s_values, allos1_q_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos2_s_values, allos2_q_values, 'Color', [0 0.4470 0.7410 .6])
        plot(autos1_s_values, autos1_q_values, 'Color', [0.8500 0.3250 0.0980])
        plot(autos2_s_values, autos2_q_values, 'Color', [0.8500 0.3250 0.0980 .6])
        plot(diploids_s_values, diploids_q_values, 'Color', 'k')

        if h_val == 0
            legend('0 HEs', '2 HEs', 'α=0', 'α=1/4', 'Diploids')
        end

        if h_val == 0
            title('h=0 (recessive)')
        elseif h_val == .5
            title('h=.5 (additive)')
        elseif h_val == 1
            title('h=1 (dominant)')
        elseif h_val > .5 && h_val < 1 
            title('partial dominance')
        elseif h_val > 0 && h_val < .5
            title('partial recessivity')
        end

        xscale log
        ylim([0, 1])
        if grid_pos == 1 
            ylabel('q (deleterious allele)')
        end

        % eigenvalue plot
        subplot(grid_row, grid_col, grid_pos+3)

        plot(allos1_s_values, allos1_eig_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos2_s_values, allos2_eig_values, 'Color', [0 0.4470 0.7410 .6])
        plot(autos1_s_values, autos1_eig_values, 'Color', [0.8500 0.3250 0.0980])
        plot(autos2_s_values, autos2_eig_values, 'Color', [0.8500 0.3250 0.0980 .6])
        
        xscale log
        yscale log
        if grid_pos == 1 
            ylabel('min(|λ|)')
        end

        % stiffness ratio plot
        subplot(grid_row, grid_col, grid_pos+6)

        plot(allos1_s_values, allos1_stiff_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos2_s_values, allos2_stiff_values, 'Color', [0 0.4470 0.7410 .6])
        plot(autos1_s_values, autos1_stiff_values, 'Color', [0.8500 0.3250 0.0980])
        plot(autos2_s_values, autos2_stiff_values, 'Color', [0.8500 0.3250 0.0980 .6])

        xscale log
        yscale log
        if grid_pos == 1 
            ylabel('stiffness ratio')
        end

        
    else
        % bifn plot
        subplot(grid_row, grid_col, grid_pos)

        plot(allos1_s1, allos1_q1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos1_s2, allos1_q2, 'Color', [0 0.4470 0.7410])
        plot(allos1_s3, allos1_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')
        
        plot(allos2_s1, allos2_q1, 'Color', [0 0.4470 0.7410 .6])
        plot(allos2_s2, allos2_q2, 'Color', [0 0.4470 0.7410 .6])
        plot(allos2_s3, allos2_q3, 'Color', [0 0.4470 0.7410 .6], 'LineStyle','--')

        plot(autos1_s1, autos1_q1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s2, autos1_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s3, autos1_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(autos2_s1, autos2_q1, 'Color', [0.8500 0.3250 0.0980 .6])
        plot(autos2_s2, autos2_q2, 'Color', [0.8500 0.3250 0.0980 .6])
        plot(autos2_s3, autos2_q3, 'Color', [0.8500 0.3250 0.0980 .6], 'LineStyle','--')

        plot(diploids_s1, diploids_q1, 'Color','k')
        plot(diploids_s2, diploids_q2, 'Color','k')
        plot(diploids_s3, diploids_q3, 'Color','k', 'LineStyle','--')

        xscale log
        ylim([0, 1])
        if grid_pos == 1 
            ylabel('q (deleterious allele)')
        end

        if h_val == 0
            title('h=0 (recessive)')
        elseif h_val == .5
            title('h=.5 (additive)')
        elseif h_val == 1
            title('h=1 (dominant)')
        elseif h_val > .5 && h_val < 1 
            title('partial dominance')
        elseif h_val > 0 && h_val < .5
            title('partial recessivity')
        end

        % eigenvalue plot
        subplot(grid_row, grid_col, grid_pos+3)

        plot(allos1_s1, allos1_eig1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos1_s2, allos1_eig2, 'Color', [0 0.4470 0.7410])
       
        plot(allos2_s1, allos2_eig1, 'Color', [0 0.4470 0.7410 .6])
        plot(allos2_s2, allos2_eig2, 'Color', [0 0.4470 0.7410 .6])

        plot(autos1_s1, autos1_eig1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s2, autos1_eig2, 'Color', [0.8500 0.3250 0.0980])

        plot(autos2_s1, autos2_eig1, 'Color', [0.8500 0.3250 0.0980 .6])
        plot(autos2_s2, autos2_eig2, 'Color', [0.8500 0.3250 0.0980 .6])

        xscale log
        yscale log
        if grid_pos == 1 
            ylabel('min(|λ|)')
        end

        % stiffness ratio plot
        subplot(grid_row, grid_col, grid_pos+6)

        plot(allos1_s1, allos1_stiff1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(allos1_s2, allos1_stiff2, 'Color', [0 0.4470 0.7410])
       
        plot(allos2_s1, allos2_stiff1, 'Color', [0 0.4470 0.7410 .6])
        plot(allos2_s2, allos2_stiff2, 'Color', [0 0.4470 0.7410 .6])

        plot(autos1_s1, autos1_stiff1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s2, autos1_stiff2, 'Color', [0.8500 0.3250 0.0980])

        plot(autos2_s1, autos2_stiff1, 'Color', [0.8500 0.3250 0.0980 .6])
        plot(autos2_s2, autos2_stiff2, 'Color', [0.8500 0.3250 0.0980 .6])

        xscale log
        yscale log
        if grid_pos == 1 
            ylabel('stiffness ratio')
        end

    end

    if grid_pos+6 == 8
        xlabel('s (selection coefficient)')
    end
    
    

end
