% for autos, creates a 2D plot of q values over variable s

iterations = 100; % number of steps or number of data points to generate

s_init_val = 1.7e-4; % starting s value
s_step_size = 1e-7; % size of change in s for each iteration

mu_val = 1e-6; % constant value of mutation rate
nu_val = 1e-7; % constant value of backward mutation rate
%mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
a_val = 0; % constant value of alpha (double reduction rate)

h1_val = 1; % h1 dominance coefficient value, constant
h2_val = 1; % h2 dominance coefficient value, constant
h3_val = 1; % h3 dominance coefficient value, constant


syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
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
sel_g0 = G0*w0+(1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3;
sel_g1 = (1/2 - a/2)*G1*w1 + (2/3 - 2*a/3)*G2*w2 + (1/2 - a/2)*G3*w3;
sel_g2 = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4;

% equations for mutation
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0 == 0;
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu - g1 == 0;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2 == 0;

mut_eqn_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_eqn_set)
    mut_eqn_set(i) = subs(mut_eqn_set(i), G0, g0^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G1, 2*g0*g1);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G2, (2*g0*g2 + g1^2));
    mut_eqn_set(i) = subs(mut_eqn_set(i), G3, 2*g1*g2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G4, g2^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), g2, (1-g1-g0));
end

%creates arrays/matrices to store data
g2_values_array = ones(1, iterations);
g0_max_values_array = zeros(1, iterations);
g1_max_values_array = zeros(1, iterations);
g0_min_values_array = zeros(1, iterations);
g1_min_values_array = zeros(1, iterations);
s_values_array = zeros(1, iterations);
s_current_val = s_init_val;

%iterates through values of s to find the stable fixed point
for i = 1:iterations

    %stores the s value in an array
    s_values_array(i) = s_current_val;
    
    %calls a numeric solving function to find all fixed points
    [g0_value, g1_value] = numeric_solver(mut_eqn_set(1), mut_eqn_set(2), mu, mu_val, nu, nu_val, s, s_current_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);

    %selects the stable fixed point as being that with the largest g0
    %value (this has not been formally proven, but has support from 
    %both biological intuition and linear stability analysis
    for j = 1:length(g0_value)
        if g0_value(j) == max(g0_value)
            g0_max_values_array(i) = g0_value(j);
            g1_max_values_array(i) = g1_value(j);
        end
    end

    for j = 1:length(g0_value)
        if g0_value(j) == min(g0_value)
            g0_min_values_array(i) = g0_value(j);
            g1_min_values_array(i) = g1_value(j);
        end
    end
    
    %increase s by the specified step size
    s_current_val = s_current_val + s_step_size;

end

%calculates the value of q using g0 and g1 values
q_max_values_array = g0_max_values_array + (1/2)*g1_max_values_array;

q_min_values_array = g0_min_values_array + (1/2)*g1_min_values_array;

%calculates the g2_values by subtracting the g0 and g1 values from an array
%of ones
g2_values_array = g2_values_array - g0_max_values_array - g1_max_values_array;

%creates strings of the input parameters to be put on the graphs
iterations_str = strcat('# steps: ', string(iterations));
s_init_str = strcat('initial s: ', string(s_init_val));
s_step_size_str = strcat('s step-size: ',string(s_step_size));
h1_str = strcat('h1: ',string(h1_val));
h2_str = strcat('h2: ',string(h2_val));
h3_str = strcat('h3: ',string(h3_val));
mu_str = strcat('mu: ',string(mu_val));
a_str = strcat('alpha: ',string(a_val));

parameters_str = {'Parameters:', s_init_str, s_step_size_str, iterations_str, mu_str, h1_str, h2_str, h3_str, a_str};
dim = [0.5 0.5 0.3 0.3];

%plots a 2D figure of the stable q values over s
figure

scatter(s_values_array, q_max_values_array)
xscale log
title('Autos: MAX Allele Frequency vs. Selection Coefficient')
ylabel('q (ancestral allele frequency)')
xlabel('s (selection coefficient)')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')

figure

scatter(s_values_array, q_min_values_array)
xscale log
title('Autos: MIN Allele Frequency vs. Selection Coefficient')
ylabel('q (ancestral allele frequency)')
xlabel('s (selection coefficient)')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')


%plots the g0, g1, and g2 stable equilibria values over s
% figure
% title('Gamete Equilibria vs. s')
% annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')
% 
% subplot(1, 3, 1)
% scatter(s_values_array, g0_values_array)
% xscale log
% title('g0 vs. s')
% xlabel('s')
% ylabel('g0')
% 
% subplot(1, 3, 2)
% scatter(s_values_array, g1_values_array)
% xscale log
% title('g1 vs. s')
% xlabel('s')
% ylabel('g1')
% 
% subplot(1, 3, 3)
% scatter(s_values_array, g2_values_array)
% xscale log
% title('g2 vs. s')
% xlabel('s')
% ylabel('g2')


%function which uses vpasolve to evaluate the fixed points of the system
function [g0_value, g1_value] = numeric_solver(mut_g0_eqn, mut_g1_eqn, mu, mu_value, nu, nu_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1)

    g0_eqn = subs(mut_g0_eqn, mu, mu_value);
    g0_eqn = subs(g0_eqn, nu, nu_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mu_value);
    g1_eqn = subs(g1_eqn, nu, nu_value);
    g1_eqn = subs(g1_eqn, s, sel_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);


    [g0_value, g1_value] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);

end
