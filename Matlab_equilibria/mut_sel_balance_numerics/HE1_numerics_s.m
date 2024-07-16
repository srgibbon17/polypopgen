% for 1 HE allos, creates a 2D plot of q values over variable s

iterations = 100; % number of steps or number of data points to generate

s_init_val = 1e-7; % starting s value
s_step_size = 1e-7; % size of change in s for each iteration

mu_val = 1e-6; % constant value of mutation rate
h1_val = .0; % h1 dominance coefficient value, constant
h2_val = .0; % h2 dominance coefficient value, constant
h3_val = .0; % h3 dominance coefficient value, constant

syms g00 g01 g10 g11 s h1 h2 h3 mu

% assumptions on the parameters of the model; theoretical bounds
assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);

% equations to parameterize relative fitnesses
wbar = (1-2*s*(h1*(g00*g10+g00*g01)+h2*(g00*g11+g01*g10)+h3*(g01*g11+g10*g11))-s*(h2*(g01^2+g10^2)+g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g00 = g00^2*w0+(17/16)*g00*g01*w1+(17/16)*g00*g10*w1+(3/16)*g01^2*w2+(3/16)*g10^2*w2+(7/16)*g01*g10*w2+(7/16)*g00*g11*w2+(1/16)*g01*g11*w3+(1/16)*g10*g11*w3;
sel_g10 = (3/16)*g00*g01*w1+(11/16)*g00*g10*w1+(1/16)*g01^2*w2+(9/16)*g10^2*w2+(9/16)*g01*g10*w2+(9/16)*g00*g11*w2+(3/16)*g01*g11*w3+(11/16)*g10*g11*w3;
sel_g01 = (11/16)*g00*g01*w1+(3/16)*g00*g10*w1+(9/16)*g01^2*w2+(1/16)*g10^2*w2+(9/16)*g01*g10*w2+(9/16)*g00*g11*w2+(11/16)*g01*g11*w3+(3/16)*g10*g11*w3;
sel_g11 = (1/16)*g00*g01*w1+(1/16)*g00*g10*w1+(3/16)*g01^2*w2+(3/16)*g10^2*w2+(7/16)*g01*g10*w2+(7/16)*g00*g11*w2+(17/16)*g01*g11*w3+(17/16)*g10*g11*w3+g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 - g00 == 0;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu) - g01 == 0;
mut_g10 = sel_g00*mu*(1-mu) + sel_g10*(1-mu) - g10 == 0;
mut_g11 = sel_g00*mu^2 + sel_g01*mu + sel_g10*mu + sel_g11 - g11 == 0;

mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g11 from the equation by replacing it with 1-(g00+g01+g10)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g11, 1-(g00+g01+g10));

    % removes g10 from the equation using g01 = g10
    mut_eqn_set(i) = subs(mut_eqn_set(i), g10, g01);
end

%creates arrays/matrices to store data
g00_values_array = zeros(1, iterations);
g01_values_array = zeros(1, iterations);
s_values_array = zeros(1, iterations);

s_current_val = s_init_val;

%iterates through values of s to find the stable fixed point
for i = 1:iterations

    %stores the s value in an array
    s_values_array(i) = s_current_val;

    %calls a numeric solving function to find all fixed points
    [g00_value, g01_value] = numeric_solver(mut_eqn_set(1), mut_eqn_set(2), mu, mu_val, s, s_current_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g01);

    %selects the stable fixed point as being that with the largest g00
    %value (this has not been formally proven, but has support from 
    %both biological intuition and linear stability analysis
    for j = 1:length(g00_value)
        if g00_value(j) == max(g00_value)
            g00_values_array(i) = g00_value(j);
            g01_values_array(i) = g01_value(j);
        end
    end

    %increase s by the specified step size
    s_current_val = s_current_val + s_step_size;

end

%calculates the value of q using g00 and g01 values
q_values_array = g00_values_array + g01_values_array;

%creates strings of the input parameters to be put on the graphs
iterations_str = strcat('# steps: ', string(iterations));
s_init_str = strcat('initial s: ', string(s_init_val));
s_step_size_str = strcat('s step-size: ',string(s_step_size));
h1_str = strcat('h1: ',string(h1_val));
h2_str = strcat('h2: ',string(h2_val));
h3_str = strcat('h3: ',string(h3_val));
mu_str = strcat('mu: ',string(mu_val));

parameters_str = {'Parameters:', s_init_str, s_step_size_str, iterations_str, mu_str, h1_str, h2_str, h3_str};
dim = [0.5 0.5 0.3 0.3];

%plots a 2D figure of the stable q values over s
figure

scatter(s_values_array, q_values_array, 'filled')
xscale log
title('HE1: Allele Frequency vs. Selection Coefficient')
ylabel('q (ancestral allele frequency)')
xlabel('s (selection coefficient)')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')

%plots a 2D figure of the stable q values over s
figure

plot(s_values_array, q_values_array)
xscale log
title('HE1: Allele Frequency vs. Selection Coefficient')
ylabel('q (ancestral allele frequency)')
xlabel('s (selection coefficient)')
annotation('textbox', dim, 'String', parameters_str, 'FitBoxToText','on')


%function which uses vpasolve to evaluate the fixed points of the system
function [g00_value, g01_value] = numeric_solver(mut_g00_eqn, mut_g01_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, g00, g01)

    g00_eqn = subs(mut_g00_eqn, mu, mut_value);
    g00_eqn = subs(g00_eqn, s, sel_value);
    g00_eqn = subs(g00_eqn, h1, h1_value);
    g00_eqn = subs(g00_eqn, h2, h2_value);
    g00_eqn = subs(g00_eqn, h3, h3_value);

    g01_eqn = subs(mut_g01_eqn, mu, mut_value);
    g01_eqn = subs(g01_eqn, s, sel_value);
    g01_eqn = subs(g01_eqn, h1, h1_value);
    g01_eqn = subs(g01_eqn, h2, h2_value);
    g01_eqn = subs(g01_eqn, h3, h3_value);


    [g00_value, g01_value] = vpasolve([g00_eqn, g01_eqn], [g00, g01]);

end