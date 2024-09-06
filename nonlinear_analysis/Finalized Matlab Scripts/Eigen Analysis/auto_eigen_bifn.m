function [neutral_q, neutral_s, neutral_small_eigen, neutral_stiff_ratio, selected_q, selected_s, selected_small_eigen, selected_stiff_ratio, unstable_q, unstable_s] = auto_eigen_bifn(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val)
%Generates autopolyploid bifurcation plotting data
%   for autos, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       6 data arrays in the following order:
%           neutral_stable (for q and s, respectively)
%           selected_stable (for q and s, respectively)
%           unstable (for q and s, respectively)

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
g0_t_1 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2);
g1_t_1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu;
g2_t_1 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2);

delta_g0 = g0_t_1 - g0;
delta_g1 = g1_t_1 - g1;
delta_g2 = g2_t_1 - g2;

diff_eqn_set = [delta_g0, delta_g1, delta_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(diff_eqn_set)
    diff_eqn_set(i) = subs(diff_eqn_set(i), G0, g0^2);
    diff_eqn_set(i) = subs(diff_eqn_set(i), G1, 2*g0*g1);
    diff_eqn_set(i) = subs(diff_eqn_set(i), G2, (2*g0*g2 + g1^2));
    diff_eqn_set(i) = subs(diff_eqn_set(i), G3, 2*g1*g2);
    diff_eqn_set(i) = subs(diff_eqn_set(i), G4, g2^2);
    diff_eqn_set(i) = subs(diff_eqn_set(i), g2, (1-g1-g0));
end


%creates the Jacobian of the system
jac_matrix = [diff(diff_eqn_set(1), g0), diff(diff_eqn_set(1), g1); diff(diff_eqn_set(2), g0), diff(diff_eqn_set(2), g1)];

neutral_g0 = [];
neutral_g1 = [];
neutral_small_eigen = [];
neutral_stiff_ratio = [];
neutral_s = [];

selected_g0 = [];
selected_g1 = [];
selected_small_eigen = [];
selected_stiff_ratio = [];
selected_s = [];

unstable_g0 = [];
unstable_g1 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g0_root_vals, g1_root_vals] = root_solns(diff_eqn_set(1), diff_eqn_set(2), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities, stiff_ratios, small_eigens] = linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals); 

    for j = 1:length(fixed_pt_stabilities)

        if fixed_pt_stabilities(j) == 0
            unstable_g0(end+1) = g0_root_vals(j);
            unstable_g1(end+1) = g1_root_vals(j);
            unstable_s(end+1) = s_val_range(i);

        elseif fixed_pt_stabilities(j) == 1
            if g0_root_vals(j) > .3333
                selected_g0(end+1) = g0_root_vals(j);
                selected_g1(end+1) = g1_root_vals(j);
                selected_stiff_ratio(end+1) = stiff_ratios(j);
                selected_small_eigen(end+1) = small_eigens(j);
                selected_s(end+1) = s_val_range(i);
            else
                neutral_g0(end+1) = g0_root_vals(j);
                neutral_g1(end+1) = g1_root_vals(j);
                neutral_stiff_ratio(end+1) = stiff_ratios(j);
                neutral_small_eigen(end+1) = small_eigens(j);
                neutral_s(end+1) = s_val_range(i);
            end
        end
    end
end

neutral_q = ones(1, length(neutral_g0)) - neutral_g0 - .5*neutral_g1;
selected_q = ones(1, length(selected_g0)) - selected_g0 - .5*selected_g1;
unstable_q = ones(1, length(unstable_g0)) - unstable_g0 - .5*unstable_g1;


end

%%% FUNCTIONS %%%

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

function [fixed_pt_stabilities, stiff_ratios, small_eigens] = linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, det(J) and tr(J), classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g0_root_vals;
    stiff_ratios = g0_root_vals;
    small_eigens = g0_root_vals;

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

        eigenvalues = eig(jacobian_eval);

        [vectors, values] = eig(jacobian_eval)

        abs_eigens = abs(eigenvalues);

        stiff_ratios(i) = (max(abs_eigens)/min(abs_eigens));

        small_eigens(i) = min(abs_eigens);

        
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
