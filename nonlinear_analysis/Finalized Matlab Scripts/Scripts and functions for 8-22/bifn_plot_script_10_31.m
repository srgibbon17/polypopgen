% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 10; 

grid_row = 2;
grid_col = 3;

s_val_range_1 = logspace(-9, -4, iterations);
s_val_range_2 = sort(-logspace(-9, -4, iterations)); % set of selection coefficients

mu_val_1 = 2e-8;
nu_val_1 = 1e-9;

mu_val_2 = 1e-9; % forward mutation rate
nu_val_2 = 2e-8; % backward mutation rate

diploid_color = '#f04d13';
auto_color = '#66bed6';
allo_color = '#7dab5b';

% fully recessive case ----------------------------------------------------

k_val = 1; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 1;

ploidy_comparison_plot_2(s_val_range_1, mu_val_1, nu_val_1, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

grid_pos = 4;

ploidy_comparison_plot_2(s_val_range_2, mu_val_2, nu_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

% additive case -----------------------------------------------------------

k_val = .5; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 2;

ploidy_comparison_plot_2(s_val_range_1, mu_val_1, nu_val_1, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

grid_pos = 5;

ploidy_comparison_plot_2(s_val_range_2, mu_val_2, nu_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

% fully dominant case -----------------------------------------------------

k_val = 0; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 3;

ploidy_comparison_plot_2(s_val_range_1, mu_val_1, nu_val_1, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

grid_pos = 6;

ploidy_comparison_plot_2(s_val_range_2, mu_val_2, nu_val_2, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos, auto_color, allo_color, diploid_color)

savefig("6_bifn_diagram.fig")