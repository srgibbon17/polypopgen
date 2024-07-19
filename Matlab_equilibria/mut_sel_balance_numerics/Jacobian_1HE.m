% for 1 HE allos, classification of fixed points using linear stability
% analysis and the Jacobian matrix

iterations = 11; % number of steps or number of data points to generate

s_init_val = 1e-4; % starting s value

mu_val = 1e-6; % constant value of mutation rate
h1_val = 1; % h1 dominance coefficient value, constant
h2_val = 1; % h2 dominance coefficient value, constant
h3_val = 1; % h3 dominance coefficient value, constant

coord_range = .001; % sets range for vector field

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

% expressions for mutation
mut_g00 = sel_g00*(1-mu)^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu) - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g10*(1-mu) - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu + sel_g10*mu + sel_g11 - g11;

mut_exp_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_exp_set)
    % removes g11 from the equation by replacing it with 1-(g00+g01+g10)
    mut_exp_set(i) = subs(mut_exp_set(i), g11, 1-(g00+g01+g10));
    
    % removes g10 from the equation using g01 = g10
    mut_exp_set(i) = subs(mut_exp_set(i), g10, g01);
    
end


%solves for the fixed points of the system
[g00_value, g01_value] = numeric_solver(mut_exp_set(1), mut_exp_set(2), mu, mu_val, s, s_init_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g01);

%creates the Jacobian of the system
jacobian_1 = [diff(mut_exp_set(1), g00), diff(mut_exp_set(1), g01); diff(mut_exp_set(2), g00), diff(mut_exp_set(2), g01)];

%initializes equations to create vectors for the quiver plot
[g00_eqn, g01_eqn] = quiver_plot_init(mut_exp_set(1), mut_exp_set(2), mu, mu_val, s, s_init_val, h1, h1_val, h2, h2_val, h3, h3_val);

% for each fixed point, evaluates the jacobian at that point
% for the evaluated jacobian, calculates eigenvalues and vectors
% uses the determinant and trace to classify the fixed points of the system
for i = 1:length(g00_value)
    jacobian_eval = zeros(length(jacobian_1));
    %evaluating the jacobian
    for j = 1:length(jacobian_eval)
        for k = 1:length(jacobian_eval)
           jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, s, s_init_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_value(i), g01, g01_value(i)); 
        end
    end
    
    %calulating the trace and determinant of the evaluated jacobian
    trace_jac = trace(jacobian_eval);
    det_jac = det(jacobian_eval);

    %creates a string of the current point
    current_pt_str = strcat(string(g00_value(i)), ', ', string(g01_value(i)));

    %classifies the fixed point according to the trace and determinant
    if det_jac < 0
        disp(strcat(current_pt_str, " is a saddle point"))
    elseif det_jac == 0
        ddisp(strcat(current_pt_str, " is non-isolated"))
    elseif det_jac > 0
        if trace_jac == 0
            disp(strcat(current_pt_str, " is a center"))
        elseif trace_jac^2 - 4*det_jac == 0
            disp(strcat(current_pt_str, " is a star or degenerate node"))
        elseif trace_jac > 0 && trace_jac^2 - 4*det_jac < 0
            disp(strcat(current_pt_str, " is an unstable spiral"))
        elseif trace_jac > 0 && trace_jac^2 - 4*det_jac > 0
            disp(strcat(current_pt_str, " is an unstable node"))
        elseif trace_jac < 0 && trace_jac^2 - 4*det_jac < 0
            disp(strcat(current_pt_str, " is a stable spiral"))
        elseif trace_jac < 0 && trace_jac^2 - 4*det_jac > 0
            disp(strcat(current_pt_str, " is a stable node"))
        end
    end

    %computes the eigenvectors and values of the jacobian
    [eigenvectors, eigenvalues] = eig(jacobian_eval)

    %g00_input_values = 0:1/(iterations-1):1;
    %g01_input_values = 0:1/(iterations-1):1;

    %creates a small range of values near the fixed point
    g00_input_values = (g00_value(i)-coord_range):2*coord_range/(iterations-1):(g00_value(i)+coord_range);
    g01_input_values = (g01_value(i)-coord_range):2*coord_range/(iterations-1):(g01_value(i)+coord_range);

    %creates coordinate data for the quiver plot using meshgrid
    [g00_indexing_values, g01_indexing_values] = meshgrid(g00_input_values, g01_input_values);

    %creates a blank array to store vector values for the quiver plot
    g00_vector_values = zeros(iterations, iterations);
    g01_vector_values = zeros(iterations, iterations);

    %generates the vectors for the quiver plot
    for j = 1:iterations
        for k = 1:iterations
            [g00_vector, g01_vector] = quiver_plot_vectors(g00_eqn, g01_eqn, g00_indexing_values(j, k), g01_indexing_values(j, k), g00, g01);
            g00_vector_values(j, k) = g00_vector;
            g01_vector_values(j, k) = g01_vector;
        end
    end

    %creates a figure for the quiver plot and fixed point
    figure

    quiver(g00_indexing_values, g01_indexing_values, g00_vector_values, g01_vector_values)
    hold on
    plot(g00_value(i), g01_value(i), '.', 'MarkerSize', 10)

    for j = 1:length(eigenvectors)
        slope = eigenvectors(2, j)/eigenvectors(1, j);
        x = linspace(g00_value(i)-coord_range, g00_value(i)+coord_range, 50);
        y = slope*(x - g00_value(i)) + g01_value(i);
        plot(x, y);
    end

    title('1HE: Phase Space Flow Diagram')
    xlabel('g00 value')
    ylabel('g01 value')
end

figure

x_1 = linspace(0, 1, 100);
y_1 = zeros(1, 100);
y_2 = zeros(1, 100);

for i = 1:length(x_1)
    y_1_soln = vpasolve(subs(g00_eqn==0, g00, x_1(i)), g01);
    for j = 1:length(y_1_soln)    
        if y_1_soln(j) >= 0 && y_1_soln(j) <= 1
            y_1(i) = y_1_soln(j);
        end
    end
    y_2_soln = vpasolve(subs(g01_eqn==0, g00, x_1(i)), g01);
    for j = 1:length(y_2_soln)    
        if y_2_soln(j) >= 0 && y_2_soln(j) <= 1
            y_2(i) = y_2_soln(j);
        end
    end
end

plot(x_1, y_1)
hold on
plot(x_1, y_2)

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

%function which acts as a partial derivative evaluation tool 
%used to evaluate the jacobian at each entry
function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_value, s, s_value, h1, h1_value, h2, h2_value, h3, h3_value, g00, g00_sub_value, g01, g01_sub_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, s, s_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, g00, g00_sub_value);
    pd_value = subs(pd_value, g01, g01_sub_value);

end

%initializes the quiver plot by substituting in all of the parameter values
%which are constant (i.e. s, mu, h1, h2, h3)
function [g00_eqn, g01_eqn] = quiver_plot_init(mut_g00_eqn, mut_g01_eqn, mu, mut_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value)

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

end

%generates vectors for the quiver plot by substituting the current values
%of g00 and g01
function [g00_vector, g01_vector] = quiver_plot_vectors(g00_eqn, g01_eqn, g00_value, g01_value, g00, g01)

    g00_vector = subs(g00_eqn, g00, g00_value);
    g00_vector = subs(g00_vector, g01, g01_value);

    g01_vector = subs(g01_eqn, g00, g00_value);
    g01_vector = subs(g01_vector, g01, g01_value);

end