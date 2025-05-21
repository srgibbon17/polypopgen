s_val_range = logspace(-9, -3, 500);
mu_val = 2e-8;
nu_val = 1e-9;
h1_val = 1;
h2_val = 1;
h3_val = 1;

%%%%%%%

eta_val = 0;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(selected, 'allo_selected_eta_0.csv')
writematrix(neutral, 'allo_neutral_eta_0.csv')
writematrix(unstable, 'allo_unstable_eta_0.csv')

%%%%%%

eta_val = 1/4;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_4_kappa_0.csv')
writematrix(selected, 'allo_selected_eta_1_4_kappa_0.csv')
writematrix(unstable, 'allo_unstable_eta_1_4_kappa_0.csv')

kappa_val = -1/3;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_4_kappa_min.csv')
writematrix(selected, 'allo_selected_eta_1_4_kappa_min.csv')
writematrix(unstable, 'allo_unstable_eta_1_4_kappa_min.csv')

kappa_val = 1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_4_kappa_1.csv')
writematrix(selected, 'allo_selected_eta_1_4_kappa_1.csv')
writematrix(unstable, 'allo_unstable_eta_1_4_kappa_1.csv')

%%%%%%%

eta_val = 1/3;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_3_kappa_0.csv')
writematrix(selected, 'allo_selected_eta_1_3_kappa_0.csv')
writematrix(unstable, 'allo_unstable_eta_1_3_kappa_0.csv')

kappa_val = -1/2;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_3_kappa_min.csv')
writematrix(selected, 'allo_selected_eta_1_3_kappa_min.csv')
writematrix(unstable, 'allo_unstable_eta_1_3_kappa_min.csv')

kappa_val = 1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_3_kappa_1.csv')
writematrix(selected, 'allo_selected_eta_1_3_kappa_1.csv')
writematrix(unstable, 'allo_unstable_eta_1_3_kappa_1.csv')

%%%%%%%

eta_val = 1/2;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_2_kappa_0.csv')
writematrix(selected, 'allo_selected_eta_1_2_kappa_0.csv')
writematrix(unstable, 'allo_unstable_eta_1_2_kappa_0.csv')

kappa_val = -1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_2_kappa_min.csv')
writematrix(selected, 'allo_selected_eta_1_2_kappa_min.csv')
writematrix(unstable, 'allo_unstable_eta_1_2_kappa_min.csv')

kappa_val = 1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1_2_kappa_1.csv')
writematrix(selected, 'allo_selected_eta_1_2_kappa_1.csv')
writematrix(unstable, 'allo_unstable_eta_1_2_kappa_1.csv')

%%%%%%%

eta_val = 2/3;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_2_3_kappa_0.csv')
writematrix(selected, 'allo_selected_eta_2_3_kappa_0.csv')
writematrix(unstable, 'allo_unstable_eta_2_3_kappa_0.csv')

kappa_val = -1/2;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_2_3_kappa_min.csv')
writematrix(selected, 'allo_selected_eta_2_3_kappa_min.csv')
writematrix(unstable, 'allo_unstable_eta_2_3_kappa_min.csv')

kappa_val = 1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_2_3_kappa_1.csv')
writematrix(selected, 'allo_selected_eta_2_3_kappa_1.csv')
writematrix(unstable, 'allo_unstable_eta_2_3_kappa_1.csv')

%%%%%%%

eta_val = 3/4;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_3_4_kappa_0.csv')
writematrix(selected, 'allo_selected_eta_3_4_kappa_0.csv')
writematrix(unstable, 'allo_unstable_eta_3_4_kappa_0.csv')

kappa_val = -1/3;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_3_4_kappa_min.csv')
writematrix(selected, 'allo_selected_eta_3_4_kappa_min.csv')
writematrix(unstable, 'allo_unstable_eta_3_4_kappa_min.csv')

kappa_val = 1;

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_3_4_kappa_1.csv')
writematrix(selected, 'allo_selected_eta_3_4_kappa_1.csv')
writematrix(unstable, 'allo_unstable_eta_3_4_kappa_1.csv')

%%%%%%%

eta_val = 1;
kappa_val = 0;
disp(eta_val)

[neutral, selected, unstable] = new_allo_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(neutral, 'allo_neutral_eta_1.csv')
writematrix(selected, 'allo_selected_eta_1.csv')
writematrix(unstable, 'allo_unstable_eta_1.csv')
