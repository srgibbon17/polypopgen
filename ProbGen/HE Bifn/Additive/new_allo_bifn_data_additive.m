s_val_range = logspace(-9, -3, 500);
mu_val = 2e-8;
nu_val = 1e-9;
h1_val = .25;
h2_val = .5;
h3_val = .75;

%%%%%%%

beta_val = 0;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_0.csv')

%%%%%%%

beta_val = 1/4;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_0.csv')

gamma_val = -1/3;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_min.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_1.csv')

%%%%%%%

beta_val = 1/3;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_0.csv')

gamma_val = -1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_min.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_1.csv')

%%%%%%%

beta_val = 1/2;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_0.csv')

gamma_val = -1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_min.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_1.csv')

%%%%%%%

beta_val = 2/3;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_0.csv')

gamma_val = -1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_min.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_1.csv')

%%%%%%%

beta_val = 3/4;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_0.csv')

gamma_val = -1/3;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_min.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_1.csv')

%%%%%%%

beta_val = 1;
gamma_val = 0;
disp(beta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1.csv')
