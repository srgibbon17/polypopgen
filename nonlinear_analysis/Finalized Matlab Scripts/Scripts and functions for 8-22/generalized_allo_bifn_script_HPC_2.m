s_val_range = logspace(-7, log10(.999), 500);
mu_val = 1e-8;
nu_val = 1e-8;
h1_val = 0;
h2_val = 0;
h3_val = 0;

%%%%%%%%

alpha_val = 1/12;
[stable_data_matrix] = auto_bifn_data_matrix_form(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data_matrix, 'auto_data_alpha_1_12.csv')

slope_val = 0;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_0.csv')

slope_val = 1/6;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_1_6.csv')

slope_val = 1/3;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_1_3.csv')

slope_val = 1/2;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_1_2.csv')

slope_val = 2/3;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_2_3.csv')

slope_val = 5/6;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_5_6.csv')

slope_val = 1;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_1.csv')

slope_val = 6/5;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_6_5.csv')

slope_val = 3/2;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_3_2.csv')

slope_val = 2;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_2.csv')

slope_val = 3;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_12_slope_3.csv')


%%%%%%%%%%%

alpha_val = 1/6;
[stable_data_matrix] = auto_bifn_data_matrix_form(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data_matrix, 'auto_data_alpha_1_6.csv')

slope_val = 1/3;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_1_3.csv')

slope_val = 1/2;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_1_2.csv')

slope_val = 2/3;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_2_3.csv')

slope_val = 5/6;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_5_6.csv')

slope_val = 1;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_1.csv')

slope_val = 6/5;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_6_5.csv')

slope_val = 3/2;
[beta_val, gamma_val] = beta_gamma_calculation(slope_val, alpha_val);

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_6_slope_3_2.csv')

%%%%%%%%%%%

alpha_val = 1/4;
beta_val = 1;
gamma_val = 1;

[stable_data_matrix] = auto_bifn_data_matrix_form(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, alpha_val);
writematrix(stable_data_matrix, 'auto_data_alpha_1_4.csv')

[stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);
writematrix(stable_data_matrix, 'allo_data_alpha_1_4.csv')
