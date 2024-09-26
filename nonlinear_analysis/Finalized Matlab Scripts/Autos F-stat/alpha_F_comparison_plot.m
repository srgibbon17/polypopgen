mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate

h1_val = 1;
h2_val = 1;
h3_val = 1;

s_val = .95;
a_val_range = linspace(0, 1, 34);

F_values = zeros(1, length(a_val_range));

for i = 1:length(a_val_range)
    [neutral_F, neutral_s, selected_F, selected_s, unstable_F, unstable_s] = auto_F_stat_bifn_data(s_val, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_range(i));

    F_values(i) = selected_F;
    disp(i)
end

figure

plot(a_val_range, F_values);
hold on
plot(a_val_range, a_val_range, 'Color', 'k', 'LineStyle', '--')

xlabel 'alpha'
ylabel 'F'

xlim([0, 1])
ylim([0, 1])
