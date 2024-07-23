% for autos, classification of fixed points using linear stability
% analysis, the Jacobian matrix, and eigendirections

iterations = 11; % number of steps or number of data points to generate

s_val = [1e-7, 9.16e-6, 1e-5, 2e-5, 1.75e-4, 1e-3]; % starting s value
%s_val = 1e-5;

mu_val = 1e-6; % constant value of forward mutation rate
nu_val = 1e-7; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
a_val = 0; % constant value of alpha (double reduction rate)

h1_val = 1; % h1 dominance coefficient value, constant
h2_val = 1; % h2 dominance coefficient value, constant
h3_val = 1; % h3 dominance coefficient value, constant

syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=-1 & s<=1);
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
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0;
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu - g1;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2;

mut_exp_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), G0, g0^2);
    mut_exp_set(i) = subs(mut_exp_set(i), G1, 2*g0*g1);
    mut_exp_set(i) = subs(mut_exp_set(i), G2, (2*g0*g2 + g1^2));
    mut_exp_set(i) = subs(mut_exp_set(i), G3, 2*g1*g2);
    mut_exp_set(i) = subs(mut_exp_set(i), G4, g2^2);
    mut_exp_set(i) = subs(mut_exp_set(i), g2, (1-g1-g0));
end

%creates the Jacobian of the system
jacobian_1 = [diff(mut_exp_set(1), g0), diff(mut_exp_set(1), g1); diff(mut_exp_set(2), g0), diff(mut_exp_set(2), g1)];

%bifurcation analysis using det(eigenvalues)
[g0_bifn_value_1, g1_bifn_value_1, s_bifn_value_1] = bifn_numeric_solver(mut_exp_set(1), mut_exp_set(2), jacobian_1, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1, s, [.545, .385, 9.2e-6]);

[g0_bifn_value_2, g1_bifn_value_2, s_bifn_value_2] = bifn_numeric_solver(mut_exp_set(1), mut_exp_set(2), jacobian_1, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1, s, [.015, .215, 1.75e-4]);


figure

for h = 1:length(s_val)
    if h == 1
        disp('s<mu')
    elseif h == 2
        disp('s>mu')
    else
        disp('s>>mu')
    end

    %solves for the fixed points of the system
    [g0_value, g1_value] = numeric_solver(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);

    %initializes equations to create vectors for the quiver plot
    [g0_eqn, g1_eqn] = quiver_plot_init(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val);

    % for each fixed point, evaluates the jacobian at that point
    % for the evaluated jacobian, calculates eigenvalues and vectors
    % uses the determinant and trace to classify the fixed points of the system
    for i = 1:length(g0_value)
        jacobian_eval = zeros(length(jacobian_1));
        %evaluating the jacobian
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
            jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, nu, nu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_value(i), g1, g1_value(i)); 
            end
        end
    
        %calulating the trace and determinant of the evaluated jacobian
        trace_jac = trace(jacobian_eval);
        det_jac = det(jacobian_eval);

        %creates a string of the current point
        current_pt_str = strcat(string(g0_value(i)), ', ', string(g1_value(i)));

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
        [eigenvectors, eigenvalues] = eig(jacobian_eval);
    end


    %creates strings of the input parameters to be put on the graphs
    s_init_str = strcat('s: ', string(s_val(h)));
    h1_str = strcat('h1: ',string(h1_val));
    h2_str = strcat('h2: ',string(h2_val));
    h3_str = strcat('h3: ',string(h3_val));
    mu_str = strcat('mu: ',string(mu_val));
    nu_str = strcat('nu: ',string(nu_val));
    mut_ratio_str = strcat(['mut-ratio ', newline, ...
        '(mu/nu): '],string(mut_ratio_val));
    a_str = strcat('alpha: ', string(a_val));

    parameters_str = {'Parameters:', s_init_str, mu_str, nu_str, mut_ratio_str, h1_str, h2_str, h3_str, a_str};
    dim = [0.5 0.5 0.3 0.3];


    subplot(2, 3, h)

    x_1 = linspace(0, 1, 100);
    y_1 = zeros(1, 100);
    y_2 = zeros(1, 100);

    
    for i = 1:length(x_1)
        y_1_soln = vpasolve(subs(g0_eqn, g0, x_1(i)), g1);
        for j = 1:length(y_1_soln)    
            if y_1_soln(j) >= 0 && y_1_soln(j) <= 1 && imag(y_1_soln(j)) == 0 
                y_1(i) = y_1_soln(j);
            end
        end
        y_2_soln = vpasolve(subs(g1_eqn, g0, x_1(i)), g1);
        for j = 1:length(y_2_soln)    
            if y_2_soln(j) >= 0 && y_2_soln(j) <= 1 && imag(y_2_soln(j)) == 0 
                y_2(i) = y_2_soln(j);
            end
        end
    end

    null_1 = plot(x_1, y_1, 'DisplayName', 'g0 nullcline');
    null_1.Color = "#0072BD";
    hold on
    null_2 = plot(x_1, y_2, 'DisplayName', 'g1 nullcline');
    null_2.Color = "#D95319";

    for i = 1:length(g0_value)
        
        jacobian_eval = zeros(length(jacobian_1));
        %evaluating the jacobian
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
            jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, nu, nu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_value(i), g1, g1_value(i)); 
            end
        end

        det_jac = det(jacobian_eval);
        if det_jac < 0
            plot(g0_value(i), g1_value(i), 'square', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'DisplayName', 'saddle point')
        else
            plot(g0_value(i), g1_value(i), 'o', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'stable node')
        end

        [eigenvectors, eigenvalues] = eig(jacobian_eval);
        for j = 1:length(eigenvectors)
            slope = eigenvectors(2, j)/eigenvectors(1, j);
            x = linspace(g0_value(i)-.2, g0_value(i)+.2, 50);
            y = slope*(x - g0_value(i)) + g1_value(i);
            if det_jac < 0 && sum(eigenvalues(:, j)) < 0 
               plot(x, y, 'k', 'LineStyle', '--', 'DisplayName', 'saddle pt. stable manifold')
            else 
               plot(x, y, 'k', 'DisplayName', 'eigenvectors');
            end
        end

        %creates a small range of values near the fixed point
        g0_input_values = 0:1/(iterations-1):1;
        g1_input_values = 0:1/(iterations-1):1;

        %creates coordinate data for the quiver plot using meshgrid
        [g0_coordinates, g1_coordinates] = meshgrid(g0_input_values, g1_input_values);

        %creates a blank array to store vector values for the quiver plot
        g0_vector_values = zeros(iterations, iterations);
        g1_vector_values = zeros(iterations, iterations);

        %generates the vectors for the quiver plot
        for j = 1:iterations
            for k = 1:iterations
                [g0_vector, g1_vector] = quiver_plot_vectors(g0_eqn, g1_eqn, g0_coordinates(j, k), g1_coordinates(j, k), g0, g1);
                g0_vector_values(j, k) = g0_vector;
                g1_vector_values(j, k) = g1_vector;
            end
        end

        %lineobj_1 = streamline(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, .3, 0, [0.01, 1000000]);
        %lineobj_1.Color = "#EDB120";
        
        % lineobj_2 = streamline(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, 0, .3, [0.01, 1000000]);
        % lineobj_2.Color = "#7E2F8E";
        % 
        % lineobj_3 = streamline(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, .27, .15, [0.01, 1000000]);
        % lineobj_3.Color = "#77AC30";
        % 
        % lineobj_4 = streamline(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, .35, .4, [0.01, 1000000]);
        % lineobj_4.Color = "#4DBEEE";
        % 
        % lineobj_5 = streamline(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, .9, .05, [0.01, 1000000]);
        % lineobj_5.Color = "#A2142F";

        quiver(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, 'Color', [.7 .7 .7], 'DisplayName', 'Vector Field')

    end
    
    annotation('textbox', dim, 'String', parameters_str,'FontSize', 8, 'FitBoxToText','on')
    if h1_val == 1 && h2_val == 1 && h3_val == 1
        sgtitle('Autos: Nullclines for Dominant Allele', 'FontSize', 14)
    elseif h1_val == .25 && h2_val == .5 && h3_val == .75
        sgtitle('Autos: Nullclines for Additive Allele', 'FontSize', 14)
    elseif h1_val == 0 && h2_val == 0 && h3_val == 0
        sgtitle('Autos: Nullclines for Recessive Allele', 'FontSize', 14)
    else
        sgtitle('Autos: Nullclines for Non-Additive Allele', 'FontSize', 14)
    end

    if h == 1
        title('s < s1 (small s)', 'FontSize', 12)
    elseif h == 2
        title('s ~ s1', 'FontSize', 12)
    elseif h == 3
        title('s1 < s < s2', 'FontSize', 12)
    elseif h == 4
        title('s1 < s < s2', 'FontSize', 12)
    elseif h == 5
        title('s ~ s2', 'FontSize', 12)
    elseif h == 6
        title('s>s2 (large s)', 'FontSize', 12)
    end
    
    xlabel('g0', 'FontSize', 12)
    if h == 1
        ylabel(['Purifying Selection', newline, ...
            'g1'], 'FontSize', 12)
    end

    % if h == 4
    %     ylabel(['Positive Selection', newline, ...
    %         'g1'], 'FontSize', 14)
    % end

    xlim([0 1])
    ylim([0 1])
end

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

%function which acts as a partial derivative evaluation tool 
%used to evaluate the jacobian at each entry
function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_value, nu, nu_value, s, s_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g0_sub_value, g1, g1_sub_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, nu, nu_value);
    pd_value = subs(pd_value, s, s_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, a, a_value);
    pd_value = subs(pd_value, g0, g0_sub_value);
    pd_value = subs(pd_value, g1, g1_sub_value);

end

%initializes the quiver plot by substituting in all of the parameter values
%which are constant (i.e. s, mu, h1, h2, h3, a)
function [g0_eqn, g1_eqn] = quiver_plot_init(mut_g0_eqn, mut_g1_eqn, mu, mu_value, nu, nu_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value)

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

end

%generates vectors for the quiver plot by substituting the current values
%of g0 and g1
function [g0_vector, g1_vector] = quiver_plot_vectors(g0_eqn, g1_eqn, g0_value, g1_value, g0, g1)

    g0_vector = subs(g0_eqn, g0, g0_value);
    g0_vector = subs(g0_vector, g1, g1_value);

    g1_vector = subs(g1_eqn, g0, g0_value);
    g1_vector = subs(g1_vector, g1, g1_value);

end

%for bifurcation analysis

function [pd_value] = bifn_pd_evaluation(jacobian_entry, mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, nu, nu_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, a, a_value);

end

function [g0_bifn_value, g1_bifn_value, s_bifn_value] = bifn_numeric_solver(mut_g0_eqn, mut_g1_eqn, jacobian, mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1, s, initial_conditions)

    g0_eqn = subs(mut_g0_eqn, mu, mu_value);
    g0_eqn = subs(g0_eqn, nu, nu_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mu_value);
    g1_eqn = subs(g1_eqn, nu, nu_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);

    jacobian_eval_bifn = jacobian;

    for j = 1:length(jacobian_eval_bifn)
        for k = 1:length(jacobian_eval_bifn)
            jacobian_eval_bifn(j, k) = bifn_pd_evaluation(jacobian(j, k), mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value); 
        end
    end

    [bifn_eigenvectors, bifn_eigenvalues] = eig(jacobian_eval_bifn);

    bifn_eig_det = det(bifn_eigenvalues) == 0;
    
    [g0_bifn_value, g1_bifn_value, s_bifn_value] = vpasolve([g0_eqn, g1_eqn, bifn_eig_det], [g0, g1, s], initial_conditions);

end

