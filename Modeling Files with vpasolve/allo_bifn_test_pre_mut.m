s_val_range = logspace(-9, -3, 5);
mu_val = 1e-8;
nu_val = 0;
h1_val = .25;
h2_val = .5;
h3_val = .75;

%%%%%%%

beta_val = 0;
gamma_val = 0;
disp(beta_val)

[stable_data] = allo_bifn_stable_only_nu_zero_pre_mut(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);