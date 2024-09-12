a_val = 0;
mu_val = 2e-8;
nu_val = 1e-9;
s_val_range = logspace(-9, -4, 10);
k_val_range = linspace(0, 1, 10);

[neutral_stable_q, selected_stable_q, unstable_q, s_coord, k_coord] = auto_cusp_cat_data(s_val_range, mu_val, nu_val, k_val_range, a_val);
