% w_bar_hat calculations from auto_bifn_data_w_bar using allelic
% frequencies

% note, only use without bistability (i.e. for additive or recessive
% deleterious)

HWE_genotypes_selected = [(1-selected_stable_q).^4; 4*selected_stable_q.*(1-selected_stable_q).^3; 6*selected_stable_q.^2.*(1-selected_stable_q).^2; 4*selected_stable_q.^3.*(1-selected_stable_q); selected_stable_q.^4];
HWE_genotypes_neutral = [(1-neutral_stable_q).^4; 4*neutral_stable_q.*(1-neutral_stable_q).^3; 6*neutral_stable_q.^2.*(1-neutral_stable_q).^2; 4*neutral_stable_q.^3.*(1-neutral_stable_q); neutral_stable_q.^4];
HWE_genotypes_unstable = [(1-unstable_q).^4; 4*unstable_q.*(1-unstable_q).^3; 6*unstable_q.^2.*(1-unstable_q).^2; 4*unstable_q.^3.*(1-unstable_q); unstable_q.^4];

absolute_fitnesses_selected = [ones(1, length(selected_stable_s)); 1-h1_val*selected_stable_s; 1-h2_val*selected_stable_s; 1-h3_val*selected_stable_s; 1-selected_stable_s];
absolute_fitnesses_neutral = [ones(1, length(neutral_stable_s)); 1-h1_val*neutral_stable_s; 1-h2_val*neutral_stable_s; 1-h3_val*neutral_stable_s; 1-neutral_stable_s];
absolute_fitnesses_unstable = [ones(1, length(unstable_s)); 1-h1_val*unstable_s; 1-h2_val*unstable_s; 1-h3_val*unstable_s; 1-unstable_s];

avg_fitness_estimate_selected = zeros(1, length(selected_stable_s));
avg_fitness_estimate_neutral = zeros(1, length(neutral_stable_s));
avg_fitness_estimate_unstable = zeros(1, length(unstable_s));

for i = 1:length(avg_fitness_estimate_selected)
    avg_fitness_estimate_selected(i) = dot(HWE_genotypes_selected(:, i), absolute_fitnesses_selected(:, i));
end

for i = 1:length(avg_fitness_estimate_neutral)
    avg_fitness_estimate_neutral(i) = dot(HWE_genotypes_neutral(:, i), absolute_fitnesses_neutral(:, i));
end

for i = 1:length(avg_fitness_estimate_unstable)
    avg_fitness_estimate_unstable(i) = dot(HWE_genotypes_unstable(:, i), absolute_fitnesses_unstable(:, i));
end


% residual plot
figure
plot(selected_stable_s, selected_avg_fitness - avg_fitness_estimate_selected, 'Color', "#0072BD")
hold on
plot(neutral_stable_s, neutral_avg_fitness - avg_fitness_estimate_neutral, 'Color', "#0072BD")
plot(unstable_s, unstable_avg_fitness - avg_fitness_estimate_unstable, 'Color', "#0072BD", 'LineStyle','--')
xscale log

% comparison plot
figure
plot(avg_fitness_estimate_selected, selected_avg_fitness, 'Color', "#0072BD")
hold on
plot(avg_fitness_estimate_neutral, neutral_avg_fitness, 'Color', "#0072BD")
plot(avg_fitness_estimate_unstable, unstable_avg_fitness, 'Color', "#0072BD", 'LineStyle', '--')
plot(linspace(min([selected_avg_fitness, unstable_avg_fitness, neutral_avg_fitness]), 1, 100), linspace(min([selected_avg_fitness, unstable_avg_fitness, neutral_avg_fitness]), 1, 100), 'Color', 'k', 'Linestyle', '--')

figure
plot(selected_stable_s, selected_avg_fitness, 'Color', "#0072BD")
hold on
plot(neutral_stable_s, neutral_avg_fitness, 'Color', "#0072BD")
plot(unstable_s, unstable_avg_fitness, 'Color', "#0072BD", 'LineStyle', '--')
xscale log

figure
plot(selected_stable_s, avg_fitness_estimate_selected, 'Color', "#0072BD")
hold on
xscale log
plot(neutral_stable_s, avg_fitness_estimate_neutral, 'Color', "#0072BD")
plot(unstable_s, avg_fitness_estimate_unstable, 'Color', "#0072BD", 'LineStyle', '--')
