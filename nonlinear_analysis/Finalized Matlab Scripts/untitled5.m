
iterations = 100;

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

h_val = 1; % diploid dominance coefficient

h1_val = 1; % simplex dominance coefficient
h2_val = 1; % duplex dominance coefficient
h3_val = 1; % triplex dominance coefficient

[autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);


[HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

