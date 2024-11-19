% marginal fitness plotting script...

color_1 = '#9ec9e2';
color_2 = '#6cb0d6';
color_3 = '#3c93c2';
color_4 = '#226e9c';
color_5 = '#0d4a70';

h_val = .5;
h1_val = .25;
h2_val = .5;
h3_val = .75;
a_val = 0;

mu_val = 2e-3;
nu_val = 1e-6;

s_val_range = logspace(-9, -2, 20);

[neutral_q_dip, neutral_s_dip, neutral_avg_fitness_dip, selected_q_dip, selected_s_dip, selected_avg_fitness_dip, unstable_q_dip, unstable_s_dip, unstable_avg_fitness_dip] = diploids_bifn_data_w_bar(s_val_range, mu_val, nu_val, h_val);
disp('Diploid Data Finished')
[neutral_stable_q, neutral_stable_s, neutral_avg_fitness, neutral_pan_diseq, selected_stable_q, selected_stable_s, selected_avg_fitness, selected_pan_diseq, unstable_q, unstable_s, unstable_avg_fitness, unstable_pan_diseq] = auto_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);
disp('Auto Data Finished')

s_stable_dip = [neutral_s_dip, selected_s_dip];
q_stable_dip = [neutral_q_dip, selected_q_dip];
fitness_stable_dip = [neutral_avg_fitness_dip, selected_avg_fitness_dip];

s_stable_auto = [neutral_stable_s, selected_stable_s];
q_stable_auto = [neutral_stable_q, selected_stable_q];
fitness_stable_auto = [neutral_avg_fitness, selected_avg_fitness];

figure 

subplot(1, 2, 1)

%diploids

plot(s_stable_dip, (1)./fitness_stable_dip, 'Color', color_1, 'DisplayName', 'G_0')
hold on
plot(s_stable_dip, (1-h_val*s_stable_dip)./fitness_stable_dip, 'Color', color_3, 'DisplayName', 'G_1')
plot(s_stable_dip, (1-s_stable_dip)./fitness_stable_dip, 'Color', color_5, 'DisplayName', 'G_2')
legend
ylabel('Marginal Fitness')
xlabel('s (Selection Coefficient)')
title('Diploids')
xscale log

subplot(1, 2, 2)

%autotetraploids

plot(s_stable_auto, (1)./fitness_stable_auto, 'Color', color_1, 'DisplayName', 'G_0')
hold on
plot(s_stable_auto, (1-h1_val*s_stable_auto)./fitness_stable_auto, 'Color', color_2, 'DisplayName', 'G_1')
plot(s_stable_auto, (1-h2_val*s_stable_auto)./fitness_stable_auto, 'Color', color_3, 'DisplayName', 'G_2')
plot(s_stable_auto, (1-h3_val*s_stable_auto)./fitness_stable_auto, 'Color', color_4, 'DisplayName', 'G_3')
plot(s_stable_auto, (1-s_stable_auto)./fitness_stable_auto, 'Color', color_5, 'DisplayName', 'G_4')
legend
xlabel('s (Selection Coefficient)')
title('Autotetraploids')
xscale log