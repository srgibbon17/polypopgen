h_val = 1;
h1_val = 1;
h2_val = 1;
h3_val = 1;
a_val = 0;

mu_val = 2e-8;
nu_val = 1e-9;

s_val_range = logspace(-9, -2, 75);


[neutral_q_dip, neutral_s_dip, neutral_avg_fitness_dip, selected_q_dip, selected_s_dip, selected_avg_fitness_dip, unstable_q_dip, unstable_s_dip, unstable_avg_fitness_dip] = diploids_bifn_data_w_bar(s_val_range, mu_val, nu_val, h_val);
disp('Diploid Data Finished')
[neutral_stable_q, neutral_stable_s, neutral_avg_fitness, neutral_pan_diseq, selected_stable_q, selected_stable_s, selected_avg_fitness, selected_pan_diseq, unstable_q, unstable_s, unstable_avg_fitness, unstable_pan_diseq] = auto_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);
disp('Auto Data Finished')
%[neutral_q_allo, neutral_s_allo, neutral_avg_fitness_allo, neutral_pan_diseq_allo, selected_q_allo, selected_s_allo, selected_avg_fitness_allo, selected_pan_diseq_allo, unstable_q_allo, unstable_s_allo, unstable_avg_fitness_allo, unstable_pan_diseq_allo] = allo_bifn_data_w_bar(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
%disp('Allo Data Finished')

s_stable_dip = [neutral_s_dip, selected_s_dip];
q_stable_dip = [neutral_q_dip, selected_q_dip];
fitness_stable_dip = [neutral_avg_fitness_dip, selected_avg_fitness_dip];

s_stable_auto = [neutral_stable_s, selected_stable_s];
q_stable_auto = [neutral_stable_q, selected_stable_q];
fitness_stable_auto = [neutral_avg_fitness, selected_avg_fitness];

% s_stable_allo = [neutral_s_allo, selected_s_allo];
% q_stable_allo = [neutral_q_allo, selected_q_allo];
% fitness_stable_allo = [neutral_avg_fitness_allo, selected_avg_fitness_allo];

% figure
% 
% subplot(2, 3, 1)
% 
% plot(s_stable_dip, q_stable_dip, 'DisplayName', 'Diploid')
% hold on
% plot(s_stable_dip, q_stable_auto, 'DisplayName', 'Auto')
% xscale log
% legend
% 
% subplot(2, 3, 2)
% 
% plot(s_stable_dip, q_stable_dip, 'DisplayName', 'Diploid')
% hold on
% plot(s_stable_dip, q_stable_auto, 'DisplayName', 'Auto')
% xscale log
% legend
% 
% subplot(2, 3, 3)
% 
% plot(s_stable_dip, fitness_stable_dip, 'DisplayName', 'Diploid')
% hold on
% plot(s_stable_dip, fitness_stable_auto, 'DisplayName', 'Auto')
% xscale log
% legend
% 
% subplot(2, 3, 4)
% 
% plot(s_stable_dip, q_stable_dip./q_stable_auto, 'DisplayName', 'q ratio')
% xscale log
% legend
% 
% subplot(2, 3, 5)
% 
% plot(s_stable_dip, q_stable_dip - q_stable_auto, 'DisplayName', 'q difference')
% xscale log
% legend
% 
% subplot(2, 3, 6)
% 
% plot(s_stable_dip, fitness_stable_dip./fitness_stable_auto, 'DisplayName', 'fitness ratio')
% 
% xscale log
% legend

% figure
% 
% subplot(1, 2, 1)
% plot(s_stable_dip, q_stable_dip, 'DisplayName', 'Diploid', 'LineWidth', 1.5)
% hold on
% plot(s_stable_dip, q_stable_auto, 'DisplayName', 'Auto', 'LineWidth', 1.5)
% %plot(s_stable_dip, q_stable_allo, 'DisplayName', 'Allo')
% xscale log
% legend
% xlabel('s (Selection Coefficient)')
% ylabel('q (Derived Allele Frequency)')
% 
% subplot(1, 2, 2)
% 
% plot(s_stable_dip, fitness_stable_dip-1, 'DisplayName', 'Diploid', 'LineWidth', 1.5)
% hold on
% plot(s_stable_dip, fitness_stable_auto-1, 'DisplayName', 'Auto', 'LineWidth', 1.5)
% %plot(s_stable_dip, fitness_stable_allo, 'DisplayName', 'Auto')
% xscale log
% legend
% xlabel('s (Selection Coefficient)')
% ylabel('Δw')

figure

subplot(1, 2, 1)
plot(selected_s_dip, selected_q_dip, 'DisplayName', 'Diploid', 'LineWidth', 1.5, 'Color', 'red')
hold on
plot(neutral_s_dip, neutral_q_dip, 'LineWidth', 1.5, 'Color', 'red')
plot(unstable_s_dip, unstable_q_dip, 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '--')
plot(neutral_stable_s, neutral_stable_q, 'DisplayName', 'Auto', 'LineWidth', 1.5, 'Color', 'blue')
plot(selected_stable_s, selected_stable_q, 'LineWidth', 1.5, 'Color', 'blue')
plot(unstable_s, unstable_q, 'LineWidth', 1.5, 'Color', 'blue', 'LineStyle', '--')
%plot(s_stable_dip, q_stable_allo, 'DisplayName', 'Allo')
xscale log
legend('Diploid', '', '', 'Auto')
xlabel('s (Selection Coefficient)')
ylabel('q (Derived Allele Frequency)')

subplot(1, 2, 2)

plot(selected_s_dip, selected_avg_fitness_dip-1, 'DisplayName', 'Diploid', 'LineWidth', 1.5, 'Color', 'red')
hold on
plot(neutral_stable_s, neutral_avg_fitness-1, 'DisplayName', 'Auto', 'LineWidth', 1.5, 'Color', 'blue')
plot(neutral_s_dip, neutral_avg_fitness_dip-1, 'LineWidth', 1.5, 'Color', 'red')
plot(selected_stable_s, selected_avg_fitness-1, 'LineWidth', 1.5, 'Color', 'blue')

yscale log
%plot(s_stable_dip, q_stable_allo, 'DisplayName', 'Allo')
xscale log
legend('Diploid', '', 'Auto')
xlabel('s (Selection Coefficient)')
ylabel('Δw')