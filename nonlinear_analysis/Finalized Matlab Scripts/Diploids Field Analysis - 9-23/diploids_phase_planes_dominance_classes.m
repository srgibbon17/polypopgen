% for diploids, a nonlinear model and analysis

s_init_val = .05;

s_val = s_init_val;
mu_val = 1e-8; % constant value of forward mutation rate
nu_val = 1e-8; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = 0; % h1 dominance coefficient value, constant

syms s q G0 G1 G2 g0 g1 h mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1 & 'Real');
assume(g1>=0 & g1<=1);
assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

%%% g0 is p and g1 is q

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1/w_bar;
w1 = (1-s*h)/w_bar;
w2 = (1-s)/w_bar;

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;

%removing g1 from the equations
delta_g1 = subs(mut_g1, g0, 1-g1);

%derivative for linear stability analysis
g1_deriv = diff(delta_g1, g1);

x_index = linspace(0, 1, 100);

%%%%%%%%%%
%Plotting%
%%%%%%%%%%

figure

subplot(2, 4, 1)

[g1_roots] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1);

[fixed_pt_stabilities] = linear_stability_analysis(g1_deriv, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1, g1_roots);

delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "k", 'LineWidth', 2, 'DisplayName', 'Combined Field')
hold on

for i = 1:length(g1_roots)
    if imag(g1_roots(i)) == 0
    if fixed_pt_stabilities(i) == 0
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Unstable Equilibria')
    elseif fixed_pt_stabilities(i) == 1
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Stable Equilibria')
    end
    end
end

plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)

title("Fully Recessive")

ylabel('Δ g_1')

xlim([0, 1])
ylim([-.8e-2, .5e-2])

%%%%%%%%%%%%%%%%
subplot(2, 4, 5)

[fitness_variance, normalized_fitness_var] = calc_omega_var(s_val, h_val, x_index);

plot(x_index, fitness_variance, "Color", 'k', "LineWidth", 2)
hold on
plot(x_index, normalized_fitness_var, "Color", 'k', 'LineStyle','--', 'LineWidth',2)

xlabel('g_1')
ylabel("Fitness Variance")

ylim([0, 8e-4])

%%%%%%%%%%%%%%%%
% subplot(3, 4, 9)
% 
% omega_bar = calc_avg_fitness(s_val, h_val, x_index);
% delta_omega = calc_fitness_change(omega_bar, s_val, h_val, x_index, delta_g1_values);
% 
% plot(x_index, delta_omega, "Color", 'k', 'LineWidth', 2)
% hold on
% plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)
% 
% xlabel("g_1")
% ylabel("Δ Fitness")

%%%%%%%%%%%%%%%%
subplot(2, 4, 2)

h_val = .5;

[g1_roots] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1);

[fixed_pt_stabilities] = linear_stability_analysis(g1_deriv, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1, g1_roots);

delta_g1_values_add = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values_add, 'Color', "k", 'LineWidth', 2, 'DisplayName', 'Combined Field')
hold on

for i = 1:length(g1_roots)
    if imag(g1_roots(i)) == 0
    if fixed_pt_stabilities(i) == 0
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Unstable Equilibria')
    elseif fixed_pt_stabilities(i) == 1
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Stable Equilibria')
    end
    end
end

plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)

xlim([0, 1])
ylim([-.8e-2, .5e-2])

title("Additive")

%%%%%%%%%%%%%%%%
subplot(2, 4, 6)

[fitness_variance, normalized_fitness_var] = calc_omega_var(s_val, h_val, x_index);

plot(x_index, fitness_variance, "Color", 'k', "LineWidth", 2)
hold on
plot(x_index, normalized_fitness_var, "Color", 'k', 'LineStyle','--', 'LineWidth',2)

xlabel('g_1')

ylim([0, 8e-4])

%%%%%%%%%%%%%%%%
% subplot(3, 4, 10)
% 
% omega_bar = calc_avg_fitness(s_val, h_val, x_index);
% delta_omega = calc_fitness_change(omega_bar, s_val, h_val, x_index, delta_g1_values);
% 
% plot(x_index, delta_omega, "Color", 'k', 'LineWidth', 2)
% hold on
% plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)
% 
% xlabel("g_1")

%%%%%%%%%%%%%%%%
subplot(2, 4, 3)

h_val = .98;

[g1_roots] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1);

[fixed_pt_stabilities] = linear_stability_analysis(g1_deriv, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1, g1_roots);

delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "k", 'LineWidth', 2, 'DisplayName', 'Combined Field')
hold on

for i = 1:length(g1_roots)
    if imag(g1_roots(i)) == 0
    if fixed_pt_stabilities(i) == 0
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Unstable Equilibria')
    elseif fixed_pt_stabilities(i) == 1
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Stable Equilibria')
    end
    end
end

plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)

xlim([0, 1])
ylim([-.8e-2, .5e-2])

title("Partially Dominant")

%%%%%%%%%%%%%%%%
subplot(2, 4, 7)

[fitness_variance, normalized_fitness_var] = calc_omega_var(s_val, h_val, x_index);

plot(x_index, fitness_variance, "Color", 'k', "LineWidth", 2)
hold on
plot(x_index, normalized_fitness_var, "Color", 'k', 'LineStyle','--', 'LineWidth',2)

xlabel('g_1')

ylim([0, 8e-4])

%%%%%%%%%%%%%%%%
% subplot(3, 4, 11)
% 
% omega_bar = calc_avg_fitness(s_val, h_val, x_index);
% delta_omega = calc_fitness_change(omega_bar, s_val, h_val, x_index, delta_g1_values);
% 
% plot(x_index, delta_omega, "Color", 'k', 'LineWidth', 2)
% hold on 
% plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)
% 
% xlabel("g_1")

%%%%%%%%%%%%%%%%
subplot(2, 4, 4)

h_val = 1;

[g1_roots] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1);

[fixed_pt_stabilities] = linear_stability_analysis(g1_deriv, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1, g1_roots);

delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "k", 'LineWidth', 2, 'DisplayName', 'Combined Field')
hold on

for i = 1:length(g1_roots)
    if imag(g1_roots(i)) == 0
    if fixed_pt_stabilities(i) == 0
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Unstable Equilibria')
    elseif fixed_pt_stabilities(i) == 1
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Stable Equilibria')
    end
    end
end

plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)

xlim([0, 1])
ylim([-.8e-2, .5e-2])

title("Fully Dominant")

%%%%%%%%%%%%%%%%
subplot(2, 4, 8)

[fitness_variance, normalized_fitness_var] = calc_omega_var(s_val, h_val, x_index);

plot(x_index, fitness_variance, "Color", 'k', "LineWidth", 2)
hold on
plot(x_index, normalized_fitness_var, "Color", 'k', 'LineStyle','--', 'LineWidth',2)

xlabel('g_1')

ylim([0, 8e-4])

%%%%%%%%%%%%%%%%
% subplot(3, 4, 12)
% 
% omega_bar = calc_avg_fitness(s_val, h_val, x_index);
% delta_omega = calc_fitness_change(omega_bar, s_val, h_val, x_index, delta_g1_values);
% 
% plot(x_index, delta_omega, "Color", 'k', 'LineWidth', 2)
% hold on
% plot(x_index, zeros(1, length(x_index)), 'Color', [.5, .5, .5], 'LineWidth', .5)
% 
% xlabel("g_1")


%%%%%%%%%%%
%Functions%
%%%%%%%%%%%

function [delta_g1_values] = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1_index, g1)

    delta_g1 = subs(delta_g1, mu, mu_val);
    delta_g1 = subs(delta_g1, nu, nu_val);
    delta_g1 = subs(delta_g1, s, s_val);
    delta_g1 = subs(delta_g1, h, h_val);

    delta_g1_values = zeros(1, length(g1_index));

    for i = 1:length(g1_index)
        delta_g1_values(i) = subs(delta_g1, g1, g1_index(i));
    end

end

function [g1_root] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1)

    delta_g1 = subs(delta_g1, mu, mu_val);
    delta_g1 = subs(delta_g1, nu, nu_val);
    delta_g1 = subs(delta_g1, s, s_val);
    delta_g1 = subs(delta_g1, h, h_val);
    
    g1_root = vpasolve(delta_g1, g1);
end

function [fixed_pt_stabilities] = linear_stability_analysis(g0_derivative, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0, g0_roots)
    
    g0_derivative = subs(g0_derivative, mu, mu_val);
    g0_derivative = subs(g0_derivative, nu, nu_val);
    g0_derivative = subs(g0_derivative, s, s_val);
    g0_derivative = subs(g0_derivative, h, h_val);
    
    fixed_pt_stabilities = g0_roots;

    for i = 1:length(fixed_pt_stabilities)
        g0_derivative_eval = subs(g0_derivative, g0, g0_roots(i));
        if g0_derivative_eval < 0
            fixed_pt_stabilities(i) = 1; %1 indicates a stable node
        elseif g0_derivative_eval > 0
            fixed_pt_stabilities(i) = 0; %0 indicates an unstable node
        else
            disp('Error. Linear stabilitiy analysis failed. Derivative is 0.')
        end
    end
end

function [fitness_variance] = calc_fitness_variance(s_val, h_val, q_vector)

G0_freq = (1-q_vector).^2;
G1_freq = 2*(1-q_vector).*q_vector;
G2_freq = q_vector.^2;

avg_fitness = 1 - (G1_freq*s_val*h_val + G2_freq*s_val);
load = 1- avg_fitness;

fitness_variance = G0_freq.*((-load).^2) + G1_freq.*((G1_freq*h_val*s_val-load).^2) + G2_freq.*((G2_freq*s_val-load).^2);

end

function [avg_fitness] = calc_avg_fitness(s_val, h_val, q_vector)

G1_freq = 2*(1-q_vector).*q_vector;
G2_freq = q_vector.^2;

avg_fitness = 1 - (G1_freq*s_val*h_val + G2_freq*s_val);

end

function [delta_fitness] = calc_fitness_change(avg_fitness, s_val, h_val, q_vector, delta_q_vector)

q_vec_t_plus_1 = q_vector + delta_q_vector;

G1_freq = 2*(1-q_vec_t_plus_1).*q_vec_t_plus_1;
G2_freq = q_vec_t_plus_1.^2;

avg_fitness_t_plus_1 = 1 - (G1_freq*s_val*h_val + G2_freq*s_val);

delta_fitness = avg_fitness_t_plus_1 - avg_fitness;

end

function [var_fitness, norm_var_fitness] = calc_omega_var(s_val, h_val, q_vector)

omega_bar = 1-2*s_val*h_val*(q_vector).*(1-q_vector)-(q_vector.^2)*s_val;

var_fitness = ((1-q_vector).^2).*(1-omega_bar).^2 + (2*(1-q_vector).*q_vector).*(1-s_val*h_val-omega_bar).^2 + (q_vector.^2).*(1-s_val-omega_bar).^2;

norm_var_fitness = var_fitness./omega_bar;

end