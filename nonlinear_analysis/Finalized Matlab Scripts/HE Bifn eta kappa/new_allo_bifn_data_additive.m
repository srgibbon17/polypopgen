s_val_range = logspace(-9, -3, 500);
mu_val = 2e-8;
nu_val = 1e-9;
h1_val = .25;
h2_val = .5;
h3_val = .75;

%%%%%%%

eta_val = 0;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_0.csv')

%%%%%%%

eta_val = 1/4;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_4_kappa_0.csv')

kappa_val = -1/3;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_4_kappa_min.csv')

kappa_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_4_kappa_1.csv')

%%%%%%%

eta_val = 1/3;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_3_kappa_0.csv')

kappa_val = -1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_3_kappa_min.csv')

kappa_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_3_kappa_1.csv')

%%%%%%%

eta_val = 1/2;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_2_kappa_0.csv')

kappa_val = -1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_2_kappa_min.csv')

kappa_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1_2_kappa_1.csv')

%%%%%%%

eta_val = 2/3;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_2_3_kappa_0.csv')

kappa_val = -1/2;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_2_3_kappa_min.csv')

kappa_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_2_3_kappa_1.csv')

%%%%%%%

eta_val = 3/4;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_3_4_kappa_0.csv')

kappa_val = -1/3;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_3_4_kappa_min.csv')

kappa_val = 1;

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_3_4_kappa_1.csv')

%%%%%%%

eta_val = 1;
kappa_val = 0;
disp(eta_val)

[stable_data] = new_allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, eta_val, kappa_val);
writematrix(stable_data, 'allo_eta_1.csv')
