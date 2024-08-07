% for autos, creates a vector field in the g0, g2 phase plane

iterations = 21; % number of steps or number of data points to generate

s_init_val = 1e-5; % starting s value
s_step_size = 1e-7; % size of change in s for each iteration

mu_val = 1e-6; % constant value of mutation rate
a_val = 0; % constant value of alpha (double reduction rate)

h1_val = .25; % h1 dominance coefficient value, constant
h2_val = .5; % h2 dominance coefficient value, constant
h3_val = .75; % h3 dominance coefficient value, constant

point = [.5, .5];

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
mut_g0 = sel_meiosis_g0*(1-mu)^2 - g0;
mut_g1 = 2*sel_meiosis_g0*(1-mu)*mu + sel_meiosis_g1*(1-mu) - g1;
mut_g2 = sel_meiosis_g0*mu^2 + sel_meiosis_g1*mu + sel_meiosis_g2 - g2;

mut_eqn_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_eqn_set)
    mut_eqn_set(i) = subs(mut_eqn_set(i), G0, g0^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G1, 2*g0*g1);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G2, (2*g0*g2 + g1^2));
    mut_eqn_set(i) = subs(mut_eqn_set(i), G3, 2*g1*g2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), G4, g2^2);
    mut_eqn_set(i) = subs(mut_eqn_set(i), g1, (1-g2-g0));
end

s_current_val = s_init_val;


%initializes equations which are then used to calculate the vectors for the
%quiver plot
[g0_exp, g2_exp] = quiver_plot_init(mut_eqn_set(1), mut_eqn_set(3), mu, mu_val, s, s_current_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val);

%finds all fixed points of the system
[g0_value, g2_value] = numeric_solver(mut_eqn_set(1), mut_eqn_set(3), mu, mu_val, s, s_current_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g2);

%selects the stable fixed point as being that with the largest g00
%value (this has not been formally proven, but has support from 
%both biological intuition and linear stability analysis
for j = 1:length(g0_value)
    if g0_value(j) > 0 && g0_value(j)<=1
        g0_value_eq= g0_value(j);
        g2_value_eq = g2_value(j);
    end
end


%%% sets the range of input values for g0 and g1
%essentially defines the coordinate range for the vector field


%g0_input_values = (g0_value_eq-.0001):(.0002/(iterations-1)):(g0_value_eq+.0001);
%g2_input_values = (g2_value_eq-.0001):(.0002/(iterations-1)):(g2_value_eq+.0001);

g0_input_values = 0:1/(iterations-1):1;
g2_input_values = 0:1/(iterations-1):1;

%g0_input_values = 0:.0001/(iterations-1):.0001;
%g2_input_values = 0:.001/(iterations-1):.001;
%%%

%creates coordinate data for go and g1 using meshgrid
[g0_indexing_values, g2_indexing_values] = meshgrid(g0_input_values, g2_input_values);

%creates blank arrays which will store the vector values
g0_vector_values = zeros(iterations, iterations);
g2_vector_values = zeros(iterations, iterations);

%calculates vectors by substituting in values of g0 and g2
for i = 1:iterations
    for j = 1:iterations
        [g0_vector, g2_vector] = quiver_plot_vectors(g0_exp, g2_exp, g0_indexing_values(i, j), g2_indexing_values(i, j), g0, g2);
        g0_vector_values(i, j) = g0_vector;
        g2_vector_values(i, j) = g2_vector;
    end
end
%calculates a stream which is used to plot trajectories in the phase plane
stream_1 = stream2(g0_indexing_values, g2_indexing_values, g0_vector_values, g2_vector_values, point(1), point(2), [0.1, 100000000]);

%creates the quiver plot using the vector and coordinate data
figure
quiver(g0_indexing_values, g2_indexing_values, g0_vector_values, g2_vector_values)
streamline(stream_1)
xlabel('g0')
ylabel('g2')
title('Autos: Phase Space Flow Diagram')

%function which uses vpasolve to evaluate the fixed points of the system
function [g0_value, g2_value] = numeric_solver(mut_g0_eqn, mut_g2_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g2)

    g0_eqn = subs(mut_g0_eqn, mu, mut_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g2_eqn = subs(mut_g2_eqn, mu, mut_value);
    g2_eqn = subs(g2_eqn, s, sel_value);
    g2_eqn = subs(g2_eqn, h1, h1_value);
    g2_eqn = subs(g2_eqn, h2, h2_value);
    g2_eqn = subs(g2_eqn, h3, h3_value);
    g2_eqn = subs(g2_eqn, a, a_value);


    [g0_value, g2_value] = vpasolve([g0_eqn, g2_eqn], [g0, g2]);

end

%initializes the quiver plot by substituting in all of the parameter values
%which are constant (i.e. s, mu, h1, h2, h3, a)
function [g0_eqn, g2_eqn] = quiver_plot_init(mut_g0_eqn, mut_g2_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value)

    g0_eqn = subs(mut_g0_eqn, mu, mut_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g2_eqn = subs(mut_g2_eqn, mu, mut_value);
    g2_eqn = subs(g2_eqn, s, sel_value);
    g2_eqn = subs(g2_eqn, h1, h1_value);
    g2_eqn = subs(g2_eqn, h2, h2_value);
    g2_eqn = subs(g2_eqn, h3, h3_value);
    g2_eqn = subs(g2_eqn, a, a_value);

end

%generates vectors for the quiver plot by substituting the current values
%of g0 and g2
function [g0_vector, g2_vector] = quiver_plot_vectors(g0_eqn, g2_eqn, g0_value, g2_value, g0, g2)

    g0_vector = subs(g0_eqn, g0, g0_value);
    g0_vector = subs(g0_vector, g2, g2_value);

    g2_vector = subs(g2_eqn, g0, g0_value);
    g2_vector = subs(g2_vector, g2, g2_value);

end