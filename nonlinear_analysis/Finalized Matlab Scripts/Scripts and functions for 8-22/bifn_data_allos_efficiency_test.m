s_val_range = logspace(-7, -.001, 100);
mu_val = 1e-8;
nu_val = 1e-8;
h1_val = 1;
h2_val = 1;
h3_val = 1;
beta_val = 0;
gamma_val = 0;

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
