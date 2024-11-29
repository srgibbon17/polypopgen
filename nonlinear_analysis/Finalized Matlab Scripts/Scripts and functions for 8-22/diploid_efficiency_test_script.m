s_val_range = logspace(-9, -4, 1000);
mu_val = 2e-8;
nu_val = 1e-9;
h_val = 1;

[neutral_stable_q, neutral_stable_s, neutral_avg_fitness, selected_stable_q, selected_stable_s, selected_avg_fitness, unstable_q, unstable_s, unstable_avg_fitness] = diploids_bifn_data_w_bar(s_val_range, mu_val, nu_val, h_val);
