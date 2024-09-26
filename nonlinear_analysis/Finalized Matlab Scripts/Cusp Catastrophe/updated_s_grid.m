function [updated_s_val_range] = updated_s_grid(s_val_range, mu_val, nu_val, k_bistable_range, a_val)
%Generates autopolyploid bifurcation plotting data
%   for autos, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       6 data arrays in the following order:
%           neutral_stable (for q and s, respectively)
%           selected_stable (for q and s, respectively)
%           unstable (for q and s, respectively)

g0 = sym('g0');
g1 = sym('g1');
g2 = sym('g2');
G0 = sym('G0');
G1 = sym('G1');
G2 = sym('G2');
G3 = sym('G3');
G4 = sym('G4');
s = sym('s');
h1 = sym('h1');
h2 = sym('h2');
h3 = sym('h3');
mu = sym('mu');
nu = sym('nu');
a = sym('a');

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
jac_matrix = [diff(mut_exp_set(1), g0), diff(mut_exp_set(1), g1); 
              diff(mut_exp_set(2), g0), diff(mut_exp_set(2), g1)];

[s_coord, k_coord] = meshgrid(s_val_range, k_bistable_range);

%k controls the row (first matrix index) and s controls the column (second
%matrix index)

neutral_g0 = zeros(size(s_coord));
neutral_g1 = zeros(size(s_coord));

selected_g0 = zeros(size(s_coord));
selected_g1 = zeros(size(s_coord));

unstable_g0 = zeros(size(s_coord));
unstable_g1 = zeros(size(s_coord));

for h = 1:length(k_bistable_range)
    for i = 1:length(s_val_range)

        [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_bistable_range(h));

        %solves for the fixed points of the system
        [g0_root_vals, g1_root_vals] = root_solns(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);
        
        %evaluating the jacobian and stability of each fixed point
        [fixed_pt_stabilities] = linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals); 

        for j = 1:length(fixed_pt_stabilities)

            if fixed_pt_stabilities(j) == 0
                unstable_g0(h, i) = g0_root_vals(j);
                unstable_g1(h, i) = g1_root_vals(j);

            elseif fixed_pt_stabilities(j) == 1
                if g0_root_vals(j) > .3333
                    selected_g0(h, i) = g0_root_vals(j);
                    selected_g1(h, i) = g1_root_vals(j);
                    
                else
                    neutral_g0(h, i) = g0_root_vals(j);
                    neutral_g1(h, i) = g1_root_vals(j);
                   
                end
            end
        end
    end
    disp('iteration:')
    disp(h)
end

%calculating g_2 values
neutral_g2 = ones(size(neutral_g0)) - neutral_g0 - neutral_g1;
selected_g2 = ones(size(selected_g0)) - selected_g0 - selected_g1;
unstable_g2 = ones(size(unstable_g0)) - unstable_g0 - unstable_g1;

%calculating q values
neutral_q = neutral_g2 + .5*neutral_g1;
selected_q = selected_g2 + .5*selected_g1;
unstable_q = unstable_g2 + .5*unstable_g1;

updated_s_val_range = s_val_range;

% extrapolating to the bifurcation point
for h = 1:length(k_bistable_range)
    avg_unstable_q = mean(unstable_q(h, :));

    bifn_cutoff_1 = find(selected_q(h, :) == 1, 1, 'last');
    bifn_cutoff_2 = find(neutral_q(h, :) == 1, 1, 'first');

    [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_bistable_range(h));
    
    % if no bistability occurs, creates a smooth, single plot of the
    % "neutral" and "selected" stable equilibria
    if avg_unstable_q == 1
    else

        %approximating the values of q for the bifn points:
        g0_bifn_estimate_1 = (selected_g0(h, 1+bifn_cutoff_1) + unstable_g0(h, 1+bifn_cutoff_1))/2;
        g0_bifn_estimate_2 = (neutral_g0(h, -1+bifn_cutoff_2) + unstable_g0(h, -1+bifn_cutoff_2))/2;

        g1_bifn_estimate_1 = (selected_g1(h, 1+bifn_cutoff_1) + unstable_g1(h, 1+bifn_cutoff_1))/2;
        g1_bifn_estimate_2 = (neutral_g1(h, -1+bifn_cutoff_2) + unstable_g1(h, -1+bifn_cutoff_2))/2;

        s_bifn_estimate_1 = s_coord(h, 1+bifn_cutoff_1);
        s_bifn_estimate_2 = s_coord(h, -1+bifn_cutoff_2);

        bifn_guess_1 = [g0_bifn_estimate_1, g1_bifn_estimate_1, s_bifn_estimate_1];
        bifn_guess_2 = [g0_bifn_estimate_2, g1_bifn_estimate_2, s_bifn_estimate_2];
    
        [s_bifn_value_1, residuals_1] = bifn_numeric_solver(mut_exp_set(1), mut_exp_set(2), jac_matrix, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1, s, bifn_guess_1);
        [s_bifn_value_2, residuals_2] = bifn_numeric_solver(mut_exp_set(1), mut_exp_set(2), jac_matrix, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1, s, bifn_guess_2);

        disp(residuals_1)
        disp(residuals_2)

        updated_s_val_range(end+1) = s_bifn_value_1;
        updated_s_val_range(end+1) = s_bifn_value_2;

    end
end

updated_s_val_range = sort(updated_s_val_range);

end

%%% FUNCTIONS %%%

%%% for solving the system of ODEs for fixed points/equilibria
function [g0_root_vals, g1_root_vals] = root_solns(mut_g0_eqn, mut_g1_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1)

    %function which uses vpasolve to find the fixed points/roots of the system

    g0_eqn = subs(mut_g0_eqn, mu, mu_val);
    g0_eqn = subs(g0_eqn, nu, nu_val);
    g0_eqn = subs(g0_eqn, s, s_val);
    g0_eqn = subs(g0_eqn, h1, h1_val);
    g0_eqn = subs(g0_eqn, h2, h2_val);
    g0_eqn = subs(g0_eqn, h3, h3_val);
    g0_eqn = subs(g0_eqn, a, a_val);

    g1_eqn = subs(mut_g1_eqn, mu, mu_val);
    g1_eqn = subs(g1_eqn, nu, nu_val);
    g1_eqn = subs(g1_eqn, s, s_val);
    g1_eqn = subs(g1_eqn, h1, h1_val);
    g1_eqn = subs(g1_eqn, h2, h2_val);
    g1_eqn = subs(g1_eqn, h3, h3_val);
    g1_eqn = subs(g1_eqn, a, a_val);


    [g0_root_vals, g1_root_vals] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);
end

%%% evaluating the partial derivative/Jacobian entry at a point
function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_val, g1, g1_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, s, s_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, a, a_val);
    pd_value = subs(pd_value, g0, g0_root_val);
    pd_value = subs(pd_value, g1, g1_root_val);
end

%%% for determining fixed point stability with linearization and Jacobian
function [fixed_pt_stabilities] = linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, det(J) and tr(J), classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g0_root_vals;

    jacobian_eval = [0, 0; 0, 0];

    for i = 1:length(g0_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = pd_evaluation(jacobian_matrix(j, k), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals(i), g1, g1_root_vals(i)); 
            end
        end

        %calulating the trace and determinant of the evaluated jacobian
        trace_jac = trace(jacobian_eval);
        det_jac = det(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if det_jac < 0
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        elseif det_jac > 0
            if trace_jac < 0 && trace_jac^2 - 4*det_jac > 0
                fixed_pt_stabilities(i) = 1; %1 for stable node
            else
                disp('Error. Unexpected stability type from linear stability analysis.')
            end
        else
            disp('Error. Unexpected stability type from linear stability analysis.')
        end
    end
end

%%% for extracting dominance coefficients from hill eqn. parameter
function [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val)
    
    % evaluates the three dominance coefficients for various values of k
    % assuming d = 1 under the Kacser and Burns model

    h1_val = (1/3) / (1/3 + (k_val/(1-k_val)));

    h2_val = 1 / (1 + (k_val/(1-k_val)));

    h3_val = 3 / (3 + (k_val/(1-k_val)));

end


%%% for heat map/color of avg fitness on cusp-cat diagram
function [avg_fitness] = calc_avg_fitness(g0_vals, g1_vals, g2_vals, s_coord, k_coord, s_val_range, k_val_range)

    avg_fitness = zeros(size(g0_vals));

    for i = 1:length(k_val_range)
        for j = 1:length(s_val_range)
            [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_coord(i, j));
            avg_fitness(i, j) = 1 - s_coord(i, j) * ( h1_val*2*g0_vals(i, j)*g1_vals(i, j) + h2_val*(2*g0_vals(i, j)*g2_vals(i, j) + g1_vals(i, j)*g1_vals(i,j)) + h3_val*2*g1_vals(i,j)*g2_vals(i,j) + g2_vals(i,j)*g2_vals(i,j));
        end
    end

end

function [pd_value] = bifn_pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, a, a_val);
end

%%% for solving for the bifn point numerically
function [s_bifn_value, error_estimates] = bifn_numeric_solver(mut_g0_eqn, mut_g1_eqn, jacobian, mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1, s, initial_guess)

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

    bifn_eigenvalues = eig(jacobian_eval_bifn);

    bifn_eig_product = bifn_eigenvalues(1)*bifn_eigenvalues(2) == 0;
    
    [g0_bifn_value, g1_bifn_value, s_bifn_value] = vpasolve([g0_eqn, g1_eqn, bifn_eig_product], [g0, g1, s], initial_guess);

    error_estimates = [g0_bifn_value, g1_bifn_value, s_bifn_value] - initial_guess;

end

