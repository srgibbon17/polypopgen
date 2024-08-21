%plotting script to generate a figure with 3x3 subplots which are each a
%bifurcation diagram

%recessive case:

iterations = 20; 

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 1e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

h_val = 0; % diploid dominance coefficient

h1_val = 0; % simplex dominance coefficient
h2_val = 0; % duplex dominance coefficient
h3_val = 0; % triplex dominance coefficient

[HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

[autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);

[diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

HE0_q_values = cat(2, HE0_q1, HE0_q2);
HE0_s_values = cat(2, HE0_s1, HE0_s2);

autos_q_values = cat(2, autos_q1, autos_q2);
autos_s_values = cat(2, autos_s1, autos_s2);

diploids_q_values = cat(2, diploids_q1, diploids_q2);
diploids_s_values = cat(2, diploids_s1, diploids_s2);

figure

plot(HE0_s_values, HE0_q_values, 'Color', "#0072BD")
hold on
%plot(autos_s_values, autos_q_values, 'Color', "#D95319")
plot(diploids_s_values, diploids_q_values, 'Color', "k")

title('h=0 (recessive)')
xscale log
ylabel('q (deleterious allele)')

