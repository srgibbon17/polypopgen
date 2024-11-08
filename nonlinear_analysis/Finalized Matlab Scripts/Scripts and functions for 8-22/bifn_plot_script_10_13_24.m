% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 10; 

grid_row = 1;
grid_col = 3;

s_val_range = -logspace(-9, -3, iterations); % set of selection coefficients

mu_val = 1e-9; % forward mutation rate
nu_val = 2e-8; % backward mutation rate
a_val = 0; % double reduction rate

% fully recessive case ----------------------------------------------------

k_val = 1; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 1;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% additive case -----------------------------------------------------------

k_val = .5; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 2;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% fully dominant case -----------------------------------------------------

k_val = 0; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 3;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

savefig("del_bifn_diagram_10_13_24.fig")