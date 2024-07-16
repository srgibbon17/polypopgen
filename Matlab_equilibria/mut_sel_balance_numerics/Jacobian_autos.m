% for autos, classification of fixed points using linear stability
% analysis and the Jacobian matrix

iterations = 10; % number of steps or number of data points to generate

s_init_val = .8; % starting s value
s_step_size = 1e-6; % size of change in s for each iteration

mu_val = 1e-6; % constant value of mutation rate
a_val = 0; % constant value of alpha (double reduction rate)

h1_val = 1; % h1 dominance coefficient value, constant
h2_val = 1; % h2 dominance coefficient value, constant
h3_val = 1; % h3 dominance coefficient value, constant


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
    mut_eqn_set(i) = subs(mut_eqn_set(i), g2, (1-g1-g0));
end

%solves for the fixed points of the system
[g0_value, g1_value] = numeric_solver(mut_eqn_set(1), mut_eqn_set(2), mu, mu_val, s, s_init_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);

%creates the Jacobian of the system
jacobian_1 = [diff(mut_eqn_set(1), g0), diff(mut_eqn_set(1), g1); diff(mut_eqn_set(2), g0), diff(mut_eqn_set(2), g1)];

% for each fixed point, evaluates the jacobian at that point
% for the evaluated jacobian, calculates eigenvalues and vectors
% uses the determinant and trace to classify the fixed points of the system
for i = 1:length(g0_value)
    jacobian_eval = zeros(length(jacobian_1));
    %evaluating the jacobian
    for j = 1:length(jacobian_eval)
        for k = 1:length(jacobian_eval)
           jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, s, s_init_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_value(i), g1, g1_value(i)); 
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
        disp(strcat(current_pt_str, " is non-isolated"))
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

    %calculates the eigenvalues and vectors for the jacobian
    [eigenvectors, eigenvalues] = eig(jacobian_eval);

end

%function which uses vpasolve to evaluate the fixed points of the system
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

%function which acts as a partial derivative evaluation tool 
%used to evaluate the jacobian at each entry
function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_value, s, s_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g0_sub_value, g1, g1_sub_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, s, s_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, a, a_value);
    pd_value = subs(pd_value, g0, g0_sub_value);
    pd_value = subs(pd_value, g1, g1_sub_value);

end