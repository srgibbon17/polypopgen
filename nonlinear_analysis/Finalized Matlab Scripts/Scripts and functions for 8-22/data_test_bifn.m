
k_val = .25;

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

iterations = 100; 

grid_row = 3;
grid_col = 5;

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

a_val_1 = 0;
a_val_2 = 1/12;
a_val_3 = 1/6;


[neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
