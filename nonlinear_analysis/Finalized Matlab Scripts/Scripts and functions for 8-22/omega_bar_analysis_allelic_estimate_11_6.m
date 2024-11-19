% w_bar_hat calculations from auto_bifn_data_w_bar using allelic
% frequencies

% note, only use without bistability (i.e. for additive or recessive
% deleterious or additive or dominant beneficial)

stable_s = [neutral_stable_s, selected_stable_s];
stable_q = [neutral_stable_q, selected_stable_q];
stable_avg_fitness = [neutral_avg_fitness, selected_avg_fitness];
stable_pan_diseq = [neutral_pan_diseq, selected_pan_diseq];

HWE_genotypes = [(1-stable_q).^4; 4*stable_q.*(1-stable_q).^3; 6*stable_q.^2.*(1-stable_q).^2; 4*stable_q.^3.*(1-stable_q); stable_q.^4];

absolute_fitnesses = [ones(1, length(stable_s)); 1-h1_val*stable_s; 1-h2_val*stable_s; 1-h3_val*stable_s; 1-stable_s];

avg_fitness_estimate = zeros(1, length(stable_s));

for i = 1:length(avg_fitness_estimate)
    avg_fitness_estimate(i) = dot(HWE_genotypes(:, i), absolute_fitnesses(:, i));
end

% residual plot
figure
plot(stable_s, stable_avg_fitness - avg_fitness_estimate)
xscale log

% comparison plot
figure
plot(avg_fitness_estimate, stable_avg_fitness)
hold on
plot(linspace(min(stable_avg_fitness), 1, 100), linspace(min(stable_avg_fitness), 1, 100), 'Color', 'k', 'Linestyle', '--')

figure
plot(stable_s, stable_avg_fitness)
xscale log

figure
plot(stable_s, avg_fitness_estimate)
xscale log

figure
plot(stable_s, stable_pan_diseq)
xscale log
