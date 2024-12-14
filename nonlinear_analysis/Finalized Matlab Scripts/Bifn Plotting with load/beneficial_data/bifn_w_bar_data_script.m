s_val_range = sort(-logspace(-9, -3, 1000));
mu_val = 1e-9;
nu_val = 2e-8;


% recessive case
h_val = 0;
h1_val = 0;
h2_val = 0;
h3_val = 0;

[dip_rec_neutral, dip_rec_selected, dip_rec_unstable] = diploids_bifn_data_w_bar_matrix(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_rec_neutral, auto_rec_selected, auto_rec_unstable] = auto_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_rec_neutral, allo_rec_selected, allo_rec_unstable] = allo_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_rec_neutral, 'dip_rec_neutral.csv')
writematrix(dip_rec_selected, 'dip_rec_selected.csv')
writematrix(dip_rec_unstable, 'dip_rec_unstable.csv')

writematrix(auto_rec_neutral, 'auto_rec_neutral.csv')
writematrix(auto_rec_selected, 'auto_rec_selected.csv')
writematrix(auto_rec_unstable, 'auto_rec_unstable.csv')

writematrix(allo_rec_neutral, 'allo_rec_neutral.csv')
writematrix(allo_rec_selected, 'allo_rec_selected.csv')
writematrix(allo_rec_unstable, 'allo_rec_unstable.csv')

disp('Recessive case complete')

% additive case
h_val = .5;
h1_val = .25;
h2_val = .5;
h3_val = .75;

[dip_add_neutral, dip_add_selected, dip_add_unstable] = diploids_bifn_data_w_bar_matrix(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_add_neutral, auto_add_selected, auto_add_unstable] = auto_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_add_neutral, allo_add_selected, allo_add_unstable] = allo_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_add_neutral, 'dip_add_neutral.csv')
writematrix(dip_add_selected, 'dip_add_selected.csv')
writematrix(dip_add_unstable, 'dip_add_unstable.csv')

writematrix(auto_add_neutral, 'auto_add_neutral.csv')
writematrix(auto_add_selected, 'auto_add_selected.csv')
writematrix(auto_add_unstable, 'auto_add_unstable.csv')

writematrix(allo_add_neutral, 'allo_add_neutral.csv')
writematrix(allo_add_selected, 'allo_add_selected.csv')
writematrix(allo_add_unstable, 'allo_add_unstable.csv')

disp('Additive case complete')

% dominant case
h_val = 1;
h1_val = 1;
h2_val = 1;
h3_val = 1;

[dip_dom_neutral, dip_dom_selected, dip_dom_unstable] = diploids_bifn_data_w_bar_matrix(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_dom_neutral, auto_dom_selected, auto_dom_unstable] = auto_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_dom_neutral, allo_dom_selected, allo_dom_unstable] = allo_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_dom_neutral, 'dip_dom_neutral.csv')
writematrix(dip_dom_selected, 'dip_dom_selected.csv')
writematrix(dip_dom_unstable, 'dip_dom_unstable.csv')

writematrix(auto_dom_neutral, 'auto_dom_neutral.csv')
writematrix(auto_dom_selected, 'auto_dom_selected.csv')
writematrix(auto_dom_unstable, 'auto_dom_unstable.csv')

writematrix(allo_dom_neutral, 'allo_dom_neutral.csv')
writematrix(allo_dom_selected, 'allo_dom_selected.csv')
writematrix(allo_dom_unstable, 'allo_dom_unstable.csv')

