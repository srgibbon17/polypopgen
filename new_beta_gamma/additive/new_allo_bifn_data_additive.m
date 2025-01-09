s_val_range = logspace(-7, log10(.999), 500);
mu_val = 1e-8;
nu_val = 1e-8;
h1_val = .25;
h2_val = .5;
h3_val = .75;

%%%%%%%

alpha_val = 0;
beta_val = 0;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_0.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_0.csv')

%%%%%%%

alpha_val = 1/16;
beta_val = 1/4;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_1_16.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_0.csv')

gamma_val = 1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_1_2.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_4_gamma_1.csv')

%%%%%%%

alpha_val = 1/12;
beta_val = 1/3;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_1_12.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_0.csv')

gamma_val = 1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_1_2.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_3_gamma_1.csv')

%%%%%%%

alpha_val = 1/8;
beta_val = 1/2;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_1_8.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_0.csv')

gamma_val = 1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_1_2.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1_2_gamma_1.csv')

%%%%%%%

alpha_val = 1/6;
beta_val = 2/3;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_1_6.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_0.csv')

gamma_val = 1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_1_2.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_2_3_gamma_1.csv')

%%%%%%%

alpha_val = 3/16;
beta_val = 3/4;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_3_16.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_0.csv')

gamma_val = 1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_1_2.csv')

gamma_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_3_4_gamma_1.csv')

%%%%%%%

alpha_val = 1/4;
beta_val = 1;
gamma_val = 0;
disp(alpha_val)

[stable_data] = auto_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data, 'auto_alpha_1_4.csv')

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data, 'allo_beta_1.csv')
