% for 2 HE allos, classification of fixed points using linear stability
% analysis, the Jacobian matrix, and eigendirections

iterations = 16; % number of vectors in each dimension of the vector field

s_val = [1e-8, 1e-5, .3]; % starting s value

mu_val = 1e-6; % constant value of mutation rate
h1_val = .25; % h1 dominance coefficient value, constant
h2_val = .5; % h2 dominance coefficient value, constant
h3_val = .75; % h3 dominance coefficient value, constant

coord_range = .2; % sets range for eigenvectors


syms g00 g01 g10 g11 s h1 h2 h3 mu

% assumptions on the parameters of the model; theoretical bounds
assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);
assume(s>=-1 & s<=1);
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
sel_g00 = g00^2*w0+(9/8)*g00*g01*w1+(9/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(1/8)*g01*g11*w3+(1/8)*g10*g11*w3;
sel_g10 = (3/8)*g00*g01*w1+(3/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(3/8)*g01*g11*w3+(3/8)*g10*g11*w3;
sel_g01 = (3/8)*g00*g01*w1+(3/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(3/8)*g01*g11*w3+(3/8)*g10*g11*w3;
sel_g11 = (1/8)*g00*g01*w1+(1/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(9/8)*g01*g11*w3+(9/8)*g10*g11*w3+g11^2*w4;

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


%creates the Jacobian of the system
jacobian_1 = [diff(mut_exp_set(1), g00), diff(mut_exp_set(1), g01); diff(mut_exp_set(2), g00), diff(mut_exp_set(2), g01)];

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
    [g00_value, g01_value] = numeric_solver(mut_exp_set(1), mut_exp_set(2), mu, mu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01);

    %initializes equations to create vectors for the quiver plot
    [g00_eqn, g01_eqn] = quiver_plot_init(mut_exp_set(1), mut_exp_set(2), mu, mu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val);

    % for each fixed point, evaluates the jacobian at that point
    % for the evaluated jacobian, calculates eigenvalues and vectors
    % uses the determinant and trace to classify the fixed points of the system
    for i = 1:length(g00_value)
        jacobian_eval = zeros(length(jacobian_1));
        %evaluating the jacobian
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
            jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_value(i), g01, g01_value(i)); 
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
        [eigenvectors, eigenvalues] = eig(jacobian_eval);
    end


    %creates strings of the input parameters to be put on the graphs
    s_init_str = strcat('s: ', string(s_val(h)));
    h1_str = strcat('h1: ',string(h1_val));
    h2_str = strcat('h2: ',string(h2_val));
    h3_str = strcat('h3: ',string(h3_val));
    mu_str = strcat('mu: ',string(mu_val));

    parameters_str = {'Parameters:', s_init_str, mu_str, h1_str, h2_str, h3_str};
    dim = [0.5 0.5 0.3 0.3];


    subplot(1, 3, h)

    x_1 = linspace(0, 1, 100);
    y_1 = zeros(1, 100);
    y_2 = zeros(1, 100);

    for i = 1:length(x_1)
        y_1_soln = vpasolve(subs(g00_eqn, g00, x_1(i)), g01);
        for j = 1:length(y_1_soln)    
            if y_1_soln(j) >= 0 && y_1_soln(j) <= 1
                y_1(i) = y_1_soln(j);
            end
        end
        y_2_soln = vpasolve(subs(g01_eqn, g00, x_1(i)), g01);
        for j = 1:length(y_2_soln)    
            if y_2_soln(j) >= 0 && y_2_soln(j) <= 1
                y_2(i) = y_2_soln(j);
            end
        end
    end

    null_1 = plot(x_1, y_1);
    null_1.Color = "#0072BD";
    hold on
    null_2 = plot(x_1, y_2);
    null_2.Color = "#D95319";

    for i = 1:length(g00_value)
        
        jacobian_eval = zeros(length(jacobian_1));
        %evaluating the jacobian
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
            jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, s, s_val(h), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_value(i), g01, g01_value(i)); 
            end
        end

        det_jac = det(jacobian_eval);
        if det_jac < 0
            plot(g00_value(i), g01_value(i), 'square', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
        else
            plot(g00_value(i), g01_value(i), 'o', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8)
        end

        [eigenvectors, eigenvalues] = eig(jacobian_eval);
        for j = 1:length(eigenvectors)
            slope = eigenvectors(2, j)/eigenvectors(1, j);
            x = linspace(g00_value(i)-.2, g00_value(i)+.2, 50);
            y = slope*(x - g00_value(i)) + g01_value(i);
            if det_jac < 0 && sum(eigenvalues(:, j)) < 0 
               plot(x, y, 'k', 'LineStyle', '--')
            else 
               plot(x, y, 'k');
            end
        end

        %creates a small range of values near the fixed point
        g00_input_values = 0:1/(iterations-1):1;
        g01_input_values = 0:.5/(iterations-1):.5;

        %creates coordinate data for the quiver plot using meshgrid
        [g00_coordinates, g01_coordinates] = meshgrid(g00_input_values, g01_input_values);

        %creates a blank array to store vector values for the quiver plot
        g00_vector_values = zeros(iterations, iterations);
        g01_vector_values = zeros(iterations, iterations);

        %generates the vectors for the quiver plot
        for j = 1:iterations
            for k = 1:iterations
                [g00_vector, g01_vector] = quiver_plot_vectors(g00_eqn, g01_eqn, g00_coordinates(j, k), g01_coordinates(j, k), g00, g01);
                g00_vector_values(j, k) = g00_vector;
                g01_vector_values(j, k) = g01_vector;
            end
        end

        % lineobj_1 = streamline(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, .3, 0, [0.01, 1000000]);
        % lineobj_1.Color = "#EDB120";
        % 
        % lineobj_2 = streamline(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, 0, .3, [0.01, 1000000]);
        % lineobj_2.Color = "#7E2F8E";
        % 
        % lineobj_3 = streamline(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, .27, .15, [0.01, 1000000]);
        % lineobj_3.Color = "#77AC30";
        % 
        % lineobj_4 = streamline(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, .35, .4, [0.01, 1000000]);
        % lineobj_4.Color = "#4DBEEE";
        % 
        % lineobj_5 = streamline(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, .9, .05, [0.01, 1000000]);
        % lineobj_5.Color = "#A2142F";

        quiver(g00_coordinates, g01_coordinates, g00_vector_values, g01_vector_values, 'Color', [.7 .7 .7])

    end
    
    annotation('textbox', dim, 'String', parameters_str,'FontSize', 12, 'FitBoxToText','on')
    sgtitle('2HEs: Nullclines for Dominant Allele', 'FontSize', 16)
    if h == 1
        title('s<mu', 'FontSize', 14)
    elseif h == 2
        title('s>mu', 'FontSize', 14)
    else
        title('s>>mu', 'FontSize', 14)
    end
    
    xlabel('g00', 'FontSize', 14)
    if h == 1
        ylabel('g01', 'FontSize', 14)
    end
    xlim([0 1])
    ylim([0 .5])
    if h == 3
        legend('g00 nullcline', 'g01 nullcline', 'stable node', '', '', '', '', '', '', '', '', 'saddle point', '', 'eigenvectors', 'FontSize', 12)
    end
end




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