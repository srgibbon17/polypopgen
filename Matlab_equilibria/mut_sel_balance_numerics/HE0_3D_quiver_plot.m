% for 0 HE allos, creates a 3D vector field in the g00, g01, g10 phase space

%NOTE: Currently there is a bug such that only an upper triangle of the
%vectors are being generated and plotted... not really sure why this is???

iterations = 11; % number of steps or number of data points to generate

s_init_val = 1e-3; % starting s value
s_step_size = 1e-7; % size of change in s for each iteration

mu_val = 1e-6; % constant value of mutation rate
h1_val = .25; % h1 dominance coefficient value, constant
h2_val = .5; % h2 dominance coefficient value, constant
h3_val = .75; % h3 dominance coefficient value, constant

%starting point for a trajectory in the phase plane
point = [.4, .3, .3];

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
sel_g00 = g00^2*w0+g00*g01*w1+g00*g10*w1+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2;
sel_g10 = g00*g10*w1+g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g10*g11*w3;
sel_g01 = g00*g01*w1+g01^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3;
sel_g11 = (1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3+g10*g11*w3+g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu) - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g10*(1-mu) - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu + sel_g10*mu + sel_g11 - g11;

mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g11 from the equation by replacing it with 1-(g00+g01+g10)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g11, 1-(g00+g01+g10));
end

s_current_val = s_init_val;

%initializes equations which are then used to calculate the vectors for the
%quiver plot
[g00_exp, g01_exp, g10_exp] = quiver_plot_init(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(3), mu, mu_val, s, s_current_val, h1, h1_val, h2, h2_val, h3, h3_val);

%%% sets the range of input values for g0 and g1
%essentially defines the coordinate range for the vector field

g00_input_values = 0:1/(iterations-1):1;
g01_input_values = 0:1/(iterations-1):1;
g10_input_values = 0:1/(iterations-1):1;

%g00_input_values = 0:.001/(iterations-1):.001;
%g01_input_values = 0:.001/(iterations-1):.001;
%g10_input_values = 0:.001/(iterations-1):.001;
%%%

%creates coordinate data for go and g1 using meshgrid
[g00_indexing_values, g01_indexing_values, g10_indexing_values] = meshgrid(g00_input_values, g01_input_values, g10_input_values);

%creates blank arrays which will store the vector values
g00_vector_values = zeros(iterations, iterations, iterations);
g01_vector_values = zeros(iterations, iterations, iterations);
g10_vector_values = zeros(iterations, iterations, iterations);

%calculates vectors by substituting in values of g00, g01, and g10
for i = 1:iterations
    for j = i:iterations
        for k = i:iterations
            [g00_vector, g01_vector, g10_vector] = quiver_plot_vectors(g00_exp, g01_exp, g10_exp, g00_indexing_values(i, j, k), g01_indexing_values(i, j, k), g10_indexing_values(i, j, k), g00, g01, g10);
            g00_vector_values(i, j, k) = g00_vector;
            g01_vector_values(i, j, k) = g01_vector;
            g10_vector_values(i, j, k) = g10_vector;
        end 
    end
end

%calculates a stream which is used to plot a trajectory in the phase space
streamline_data = stream3(g00_indexing_values, g01_indexing_values, g10_indexing_values, g00_vector_values, g01_vector_values, g10_vector_values, point(1), point(2), point(3));

%creates the quiver plot using the vector and coordinate data
figure
quiver3(g00_indexing_values, g01_indexing_values, g10_indexing_values, g00_vector_values, g01_vector_values, g10_vector_values)
streamline(streamline_data)
xlabel('g00')
ylabel('g01')
zlabel('g10')
title('0 HEs: Phase Space Flow Diagram')

%initializes the quiver plot by substituting in all of the parameter values
%which are constant (i.e. s, mu, h1, h2, h3)
function [g00_eqn, g01_eqn, g10_eqn] = quiver_plot_init(mut_g00_eqn, mut_g01_eqn, mut_g10_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value)

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

    g10_eqn = subs(mut_g10_eqn, mu, mut_value);
    g10_eqn = subs(g10_eqn, s, sel_value);
    g10_eqn = subs(g10_eqn, h1, h1_value);
    g10_eqn = subs(g10_eqn, h2, h2_value);
    g10_eqn = subs(g10_eqn, h3, h3_value);
end

%generates vectors for the quiver plot by substituting the current values
%of g00, g01, and g10
function [g00_vector, g01_vector, g10_vector] = quiver_plot_vectors(g00_eqn, g01_eqn, g10_eqn, g00_value, g01_value, g10_value, g00, g01, g10)

    g00_vector = subs(g00_eqn, g00, g00_value);
    g00_vector = subs(g00_vector, g01, g01_value);
    g00_vector = subs(g00_vector, g10, g10_value);

    g01_vector = subs(g01_eqn, g00, g00_value);
    g01_vector = subs(g01_vector, g01, g01_value);
    g01_vector = subs(g01_vector, g10, g10_value);

    g10_vector = subs(g10_eqn, g00, g00_value);
    g10_vector = subs(g10_vector, g01, g01_value);
    g10_vector = subs(g10_vector, g10, g10_value);

end