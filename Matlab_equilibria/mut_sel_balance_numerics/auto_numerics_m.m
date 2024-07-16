% for autos, numerical approximation of mut-sel balance for constant mu and variable s

iterations = 100; % number of steps or number of data points to generate

mu_init_val = 1e-7; % starting s value
mu_step_size = 1e-7; % size of change in s for each iteration

s_val = 5e-5; % constant value of mutation rate
a_val = 10/60; % constant value of alpha (double reduction rate)

h1_val = .25; % h1 dominance coefficient value, constant
h2_val = .5; % h2 dominance coefficient value, constant
h3_val = .75; % h3 dominance coefficient value, constant


syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu 

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(a>=0 & a<=1/6);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);
assume(G3>=0 & G3<=1);
assume(G4>=0 & G4<=1);

% equations to parameterize relative fitnesses
wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_meiosis_g0 = G0*w0+(1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3;
sel_meiosis_g1 = (1/2 - a/2)*G1*w1 + (2/3 - 2*a/3)*G2*w2 + (1/2 - a/2)*G3*w3;
sel_meiosis_g2 = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4;

% equations for mutation
mut_g0 = sel_meiosis_g0*(1-mu)^2 - g0 == 0;
mut_g1 = 2*sel_meiosis_g0*(1-mu)*mu + sel_meiosis_g1*(1-mu) - g1 == 0;
mut_g2 = sel_meiosis_g0*mu^2 + sel_meiosis_g1*mu + sel_meiosis_g2 - g2 == 0;

mut_eqn_set = [mut_g0, mut_g1, mut_g2];

for i = 1:length(mut_eqn_set)
    mut_eqn_set(i) = subs(mut_eqn_set(i), G0, g0^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G1, 2*g0*g1);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G2, (2*g0*g2 + g1^2));
    mut_eqn_set(i) = subs(mut_eqn_set(i), G3, 2*g1*g2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G4, g2^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), g2, (1-g1-g0));
end

g2_values_array = ones(1, iterations);
g0_values_array = zeros(1, iterations);
g1_values_array = zeros(1, iterations);
mu_values_array = zeros(1, iterations);
mu_current_val = mu_init_val;

for i = 1:iterations

    mu_values_array(i) = mu_current_val;

    [g0_value, g1_value] = numeric_solver(mut_eqn_set(1), mut_eqn_set(2), mu, mu_current_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);


    for j = 1:length(g0_value)
        if g0_value(j) == max(g0_value)
            g0_values_array(i) = g0_value(j);
            g1_values_array(i) = g1_value(j);
        end
    end

    mu_current_val = mu_current_val + mu_step_size;

end

q_values_array = g0_values_array + (1/2)*g1_values_array;

g2_values_array = g2_values_array - g0_values_array - g1_values_array;

iterations_str = strcat('# steps: ', string(iterations));
mu_init_str = strcat('initial mu: ', string(mu_init_val));
mu_step_size_str = strcat('mu step-size: ',string(mu_step_size));
h1_str = strcat('h1: ',string(h1_val));
h2_str = strcat('h2: ',string(h2_val));
h3_str = strcat('h3: ',string(h3_val));
s_str = strcat('s: ',string(s_val));
a_str = strcat('alpha: ',string(a_val));

parameters_str = {'Parameters:', mu_init_str, mu_step_size_str, iterations_str, s_str, h1_str, h2_str, h3_str, a_str};
dim = [0.5 0.5 0.3 0.3];

figure

scatter(mu_values_array, q_values_array)
xscale log
title('Autos: Allele Frequency vs. Selection Coefficient')
ylabel('q (ancestral allele frequency)')
xlabel('s (selection coefficient)')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')

figure
title('Gamete Equilibria vs. s')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')

subplot(1, 3, 1)
scatter(mu_values_array, g0_values_array)
xscale log
title('g0 vs. s')
xlabel('s')
ylabel('g0')

subplot(1, 3, 2)
scatter(mu_values_array, g1_values_array)
xscale log
title('g1 vs. s')
xlabel('s')
ylabel('g1')

subplot(1, 3, 3)
scatter(mu_values_array, g2_values_array)
xscale log
title('g2 vs. s')
xlabel('s')
ylabel('g2')


function [g0_value, g1_value] = numeric_solver(mut_g0_eqn, mut_g1_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1)

    g0_eqn = subs(mut_g0_eqn, mu, mut_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mut_value);
    g1_eqn = subs(g1_eqn, s, sel_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);


    [g0_value, g1_value] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);

end

