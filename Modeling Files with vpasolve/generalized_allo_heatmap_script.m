% generalized allo heatmap data generation script

iterations = 10; % number of beta and gamma values to sample

s_val = 1e-5;

mu_val = 1e-8; % constant value of forward mutation rate
nu_val = 1e-8; % constant value of backward mutation rate

h1_val = .25; % h1 dominance coefficient value, constant
h2_val = .5; % h2 dominance coefficient value, constant
h3_val = .75; % h3 dominance coefficient value, constant

[beta_coord, gamma_coord, q_coord, g00_coord, g10_coord, g11_coord] = generalized_allo_heatmap(s_val, mu_val, nu_val, h1_val, h2_val, h3_val, iterations);