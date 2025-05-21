s_val_range = logspace(-9, -3, 50);
mu_val = 1e-8;
nu_val = 0;
h1_val = .25;
h2_val = .5;
h3_val = .75;

%%%%%%%

beta_val = 0;
gamma_val = 0;

[stable_data, fixed_q_data, fixed_sga_data, fixed_sgb_data] = allo_bifn_stable_only_nu_zero(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val);

%data strcutures: s, q, qa, qb, g00, g01, g10, g11

figure

plot(stable_data(:, 1), stable_data(:, 2), color='k', LineStyle='-', LineWidth=1.5)
xscale log
hold on
plot(fixed_q_data(:, 1), fixed_q_data(:, 2), color='k', LineStyle='--', LineWidth=1.5)
plot(fixed_sga_data(:, 1), fixed_sga_data(:, 2), color=	'k', LineStyle='--', LineWidth=1.5)
plot(fixed_sgb_data(:, 1), fixed_sgb_data(:, 2), color=	'k', LineStyle='--', LineWidth=1.5)

xlabel('s')
ylabel('q')
title('Bifurcation Diagram')
ylim([-.1, 1.1])


figure

plot(stable_data(:, 1), stable_data(:, 3), color='k', LineStyle='-', LineWidth=1.5)
xscale log
hold on
plot(fixed_q_data(:, 1), fixed_q_data(:, 3), color='k', LineStyle='--', LineWidth=1.5)
plot(fixed_sga_data(:, 1), fixed_sga_data(:, 3), color=	"#0072BD", LineStyle='--', LineWidth=1.5)
plot(fixed_sgb_data(:, 1), fixed_sgb_data(:, 3), color=	"#D95319", LineStyle='--', LineWidth=1.5)

xlabel('s')
ylabel('q_a')
title('Bifurcation Diagram')
ylim([-.1, 1.1])


figure

plot(stable_data(:, 1), stable_data(:, 4), color='k', LineStyle='-', LineWidth=1.5)
xscale log
hold on
plot(fixed_q_data(:, 1), fixed_q_data(:, 4), color='k', LineStyle='--', LineWidth=1.5)
plot(fixed_sga_data(:, 1), fixed_sga_data(:, 4), color=	"#0072BD", LineStyle='--', LineWidth=1.5)
plot(fixed_sgb_data(:, 1), fixed_sgb_data(:, 4), color=	"#D95319", LineStyle='--', LineWidth=1.5)

xlabel('s')
ylabel('q_b')
title('Bifurcation Diagram')
ylim([-.1, 1.1])
