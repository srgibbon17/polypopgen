% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

figure

iterations = 10; 

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