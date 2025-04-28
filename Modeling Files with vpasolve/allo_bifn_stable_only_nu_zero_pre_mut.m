function [stable_data, fixed_q_data, fixed_sga_data, fixed_sgb_data] = allo_bifn_stable_only_nu_zero_pre_mut(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val)
%Generates allo bifurcation plotting data
%   classifies fixed points using linear stability analysis, the Jacobian 
%   matrix, and eigenvectors

%%% Assumes nu = 0 and additive or partially recessive dominance

g00 = sym('g00'); % gamete frequency with 0 derived alleles
g01 = sym('g01'); % gamete frequency with 1 derived allele in 'b' subgenome
g10 = sym('g10'); % gamete frequency with 1 derived allele in 'a' subgenome
g11 = sym('g11'); % gamete frequency with 2 derived alleles
s = sym('s'); % selection coefficient
h1 = sym('h1'); % dominance coefficient for genotype with 1 derived allele
h2 = sym('h2'); % dominance coefficient for genotype with 2 derived alleles
h3 = sym('h3'); % dominance coefficient for genotype with 3 derived alleles
mu = sym('mu'); % forward mutation rate
nu = sym('nu'); % backward mutation rate
beta = sym('beta'); % proportion of HEs in the population
gamma = sym('gamma'); % HE interference/synergy coefficient

% check to ensure beta and gamma are within plausible ranges
if beta_val > 1 || beta_val < 0
    disp('beta_val invalid. Beta must be between 0 and 1 (inclusive).');
    return;
elseif gamma_val > 1
    disp('gamma_val invalid. Gamma must be less than or equal to 1.');
    return;
elseif gamma_val < -(min(beta_val, 1-beta_val)/max(beta_val, 1-beta_val))
    disp('gamma_val invalid given beta_val. Gamma must be greater than or equal to - min(Beta, 1-Beta)/max(Beta, 1-Beta).');
    return;
end

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

% equations for mutation
m_g00 = g00*(1-mu)^2 + g01*(1-mu)*nu + g10*(1-mu)*nu + g11*nu^2;
m_g01 = g00*mu*(1-mu) + g01*(1-mu)*(1-nu) + g10*mu*nu + g11*(1-nu)*nu;
m_g10 = g00*mu*(1-mu) + g01*mu*nu + g10*(1-mu)*(1-nu) + g11*(1-nu)*nu;
m_g11 = g00*mu^2 + g01*mu*(1-nu) + g10*mu*(1-nu) + g11*(1-nu)^2;


% equations to parameterize relative fitnesses
wbar = (1-2*s*(h1*(m_g00*m_g10+m_g00*m_g01)+h2*(m_g00*m_g11+m_g01*m_g10)+h3*(m_g01*m_g11+m_g10*m_g11))-s*(h2*(m_g01^2+m_g10^2)+m_g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g00 = m_g00^2*w0 + (.5+(beta)/16)*2*m_g01*m_g00*w1 +   (.5+(beta)/16)*2*m_g10*m_g00*w1 +     ((2*beta + beta*(1-beta)*(1-gamma))/8)*m_g01^2*w2 +     ((2*beta + beta*(1-beta)*(1-gamma))/8)*m_g10^2*w2 + (.25-(beta*(1-beta)*(1-gamma))/16)*(2*m_g00*m_g11+2*m_g01*m_g10)*w2 +      ((beta)/16)*2*m_g01*m_g11*w3 +      ((beta)/16)*2*m_g10*m_g11*w3;
sel_g10 =             ((3*beta)/16)*2*m_g01*m_g00*w1 + (.5-(5*beta)/16)*2*m_g10*m_g00*w1 + ((beta + (1-gamma)*beta^2 + beta*gamma)/8)*m_g01^2*w2 +   (1-(6*beta + beta*(1-beta)*(1-gamma))/8)*m_g10^2*w2 + (.25+(beta*(1-beta)*(1-gamma))/16)*(2*m_g00*m_g11+2*m_g01*m_g10)*w2 +    ((3*beta)/16)*2*m_g01*m_g11*w3 + (.5-(5*beta)/16)*2*m_g10*m_g11*w3;
sel_g01 =          (.5-(5*beta)/16)*2*m_g01*m_g00*w1 +    ((3*beta)/16)*2*m_g10*m_g00*w1 +   (1-(6*beta + beta*(1-beta)*(1-gamma))/8)*m_g01^2*w2 + ((beta + (1-gamma)*beta^2 + beta*gamma)/8)*m_g10^2*w2 + (.25+(beta*(1-beta)*(1-gamma))/16)*(2*m_g00*m_g11+2*m_g01*m_g10)*w2 + (.5-(5*beta)/16)*2*m_g01*m_g11*w3 +    ((3*beta)/16)*2*m_g10*m_g11*w3;
sel_g11 =               ((beta)/16)*2*m_g01*m_g00*w1 +      ((beta)/16)*2*m_g10*m_g00*w1 +     ((2*beta + beta*(1-beta)*(1-gamma))/8)*m_g01^2*w2 +     ((2*beta + beta*(1-beta)*(1-gamma))/8)*m_g10^2*w2 + (.25-(beta*(1-beta)*(1-gamma))/16)*(2*m_g00*m_g11+2*m_g01*m_g10)*w2 +   (.5+(beta)/16)*2*m_g01*m_g11*w3 +   (.5+(beta)/16)*2*m_g10*m_g11*w3 + m_g11^2*w4;

% equations for mutation
diff_g00 = sel_g00 - g00;
diff_g01 = sel_g01 - g01;
diff_g10 = sel_g10 - g10;
diff_g11 = sel_g11 - g11;

mut_exp_set = [diff_g00, diff_g01, diff_g10, diff_g11];

for i = 1:length(mut_exp_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_exp_set(i) = subs(mut_exp_set(i), g10, 1-(g00+g01+g11));    
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_exp_set(1), g00), diff(mut_exp_set(1), g01), diff(mut_exp_set(1), g11); 
                diff(mut_exp_set(2), g00), diff(mut_exp_set(2), g01), diff(mut_exp_set(2), g11); 
                diff(mut_exp_set(4), g00), diff(mut_exp_set(4), g01), diff(mut_exp_set(4), g11)];

%creates a partially evaluated vector of ODEs for constant parameters
mut_exp_set_eval = mut_exp_set;

for i = 1:length(mut_exp_set_eval)
    mut_exp_set_eval(i) = allo_diff_eq_eval(mut_exp_set(i), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val);
end

%similarly, creates a partially evaluted Jacobian for constant parameters
jac_matrix_eval = jac_matrix;

for i = 1:length(jac_matrix_eval)
    for j = 1:length(jac_matrix_eval)
        jac_matrix_eval(i, j) = allo_pd_evaluation_init(jac_matrix_eval(i,j), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val);
    end
end

stable_g00 = zeros(1, length(s_val_range));
stable_g01 = zeros(1, length(s_val_range));
stable_g11 = zeros(1, length(s_val_range));
stable_s = s_val_range;

fixed_sga_g00 = [];
fixed_sga_g01 = [];
fixed_sga_g11 = [];
fixed_sga_s = [];

fixed_sgb_g00 = [];
fixed_sgb_g01 = [];
fixed_sgb_g11 = [];
fixed_sgb_s = [];

fixed_q_g00 = [];
fixed_q_g01 = [];
fixed_q_g11 = [];
fixed_q_s = [];

for i = 1:length(s_val_range)

    disp(s_val_range(i))

    %solves for the fixed points of the system
    [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_exp_set_eval(1), mut_exp_set_eval(2), mut_exp_set_eval(4), s, s_val_range(i), g00, g01, g11);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = allo_linear_stability_analysis(jac_matrix_eval, s, s_val_range(i), g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);


    for j = 1:length(fixed_pt_stabilities) 
        if fixed_pt_stabilities(j) == 1
            stable_g00(i) = g00_root_vals(j); 
            stable_g01(i) = g01_root_vals(j); 
            stable_g11(i) = g11_root_vals(j); 
        else
            if g11_root_vals(j) == 1
                fixed_q_g00(end+1) = g00_root_vals(j);
                fixed_q_g01(end+1) = g01_root_vals(j);
                fixed_q_g11(end+1) = g11_root_vals(j);
                fixed_q_s(end+1) = s_val_range(i);
            elseif g01_root_vals(j) == 0
                fixed_sga_g00(end+1) = g00_root_vals(j);
                fixed_sga_g01(end+1) = g01_root_vals(j);
                fixed_sga_g11(end+1) = g11_root_vals(j);
                fixed_sga_s(end+1) = s_val_range(i);
            else
                fixed_sgb_g00(end+1) = g00_root_vals(j);
                fixed_sgb_g01(end+1) = g01_root_vals(j);
                fixed_sgb_g11(end+1) = g11_root_vals(j);
                fixed_sgb_s(end+1) = s_val_range(i);
            end
        end
    end

end

stable_g10 = 1 - stable_g00 - stable_g01 - stable_g11;
stable_q = stable_g11 + (1/2)*stable_g01 + (1/2)*stable(g10);
stable_qa = stable_g11 + stable_g10;
stable_qb = stable_g11 + stable_g01;

stable_data_transposed = vertcat(stable_s, stable_q, stable_qa, stable_qb, stable_g00, stable_g01, stable_g10, stable_g11);
stable_data = stable_data_transposed.';

fixed_q_g10 = 1 - fixed_q_g00 - fixed_q_g01 - fixed_q_g11;
fixed_q_q = fixed_q_g11 + (1/2)*fixed_q_g01 + (1/2)*fixed_q(g10);
fixed_q_qa = fixed_q_g11 + fixed_q_g10;
fixed_q_qb = fixed_q_g11 + fixed_q_g01;

fixed_q_data_transposed = vertcat(fixed_q_s, fixed_q_q, fixed_q_qa, fixed_q_qb, fixed_q_g00, fixed_q_g01, fixed_q_g10, fixed_q_g11);
fixed_q_data = fixed_q_data_transposed.';

fixed_sga_g10 = 1 - fixed_sga_g00 - fixed_sga_g01 - fixed_sga_g11;
fixed_sga_q = fixed_sga_g11 + (1/2)*fixed_sga_g01 + (1/2)*fixed_sga(g10);
fixed_sga_qa = fixed_sga_g11 + fixed_sga_g10;
fixed_sga_qb = fixed_sga_g11 + fixed_sga_g01;

fixed_sga_data_transposed = vertcat(fixed_sga_s, fixed_sga_q, fixed_sga_qa, fixed_sga_qb, fixed_sga_g00, fixed_sga_g01, fixed_sga_g10, fixed_sga_g11);
fixed_sga_data = fixed_sga_data_transposed.';

fixed_sgb_g10 = 1 - fixed_sgb_g00 - fixed_sgb_g01 - fixed_sgb_g11;
fixed_sgb_q = fixed_sgb_g11 + (1/2)*fixed_sgb_g01 + (1/2)*fixed_sgb(g10);
fixed_sgb_qa = fixed_sgb_g11 + fixed_sgb_g10;
fixed_sgb_qb = fixed_sgb_g11 + fixed_sgb_g01;

fixed_sgb_data_transposed = vertcat(fixed_sgb_s, fixed_sgb_q, fixed_sgb_qa, fixed_sgb_qb, fixed_sgb_g00, fixed_sgb_g01, fixed_sgb_g10, fixed_sgb_g11);
fixed_sgb_data = fixed_sgb_data_transposed.';

%stable_G1 = 2*stable_g00.*stable_g10 + 2*stable_g00.*stable_g01;
%stable_G2 = stable_g10.^2 + 2*(stable_g00.*stable_g11 + stable_g01.*stable_g10) + stable_g01.^2;
%stable_G3 = 2*stable_g11.*stable_g10 + 2*stable_g11.*stable_g01;
%stable_genotypes = [stable_g00.^2; stable_G1; stable_G2; stable_G3; stable_g11.^2];

%stable_abs_fitness = [ones(1, length(stable_s)); 1-h1_val*stable_s; 1-h2_val*stable_s; 1-h3_val*stable_s; 1-stable_s];

%stable_w_bar = zeros(1, length(stable_s));

%for i = 1:length(stable_w_bar)
%    stable_w_bar(i) = dot(stable_genotypes(:, i), stable_abs_fitness(:, i));
%end

%stable_data_transposed = vertcat(stable_s, stable_q, stable_g00, stable_g01, stable_g10, stable_g11, stable_w_bar);
%stable_data = stable_data_transposed.';

end


%%% FUNCTIONS %%%

function [diff_eqn_eval] = allo_diff_eq_eval(diff_eqn, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val)

    % substitutes all constant parameters to create simplified differential
    % equations

    diff_eqn_eval = subs(diff_eqn, mu, mu_val);
    diff_eqn_eval = subs(diff_eqn_eval, nu, nu_val);
    diff_eqn_eval = subs(diff_eqn_eval, h1, h1_val);
    diff_eqn_eval = subs(diff_eqn_eval, h2, h2_val);
    diff_eqn_eval = subs(diff_eqn_eval, h3, h3_val);
    diff_eqn_eval = subs(diff_eqn_eval, beta, beta_val);
    diff_eqn_eval = subs(diff_eqn_eval, gamma, gamma_val);

end

function [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, s, s_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, s, s_val);

    g01_eqn = subs(mut_g01_eqn, s, s_val);

    g11_eqn = subs(mut_g11_eqn, s, s_val);

    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end


function [pd_value] = allo_pd_evaluation_init(jacobian_entry, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, beta, beta_val);
    pd_value = subs(pd_value, gamma, gamma_val);

end

function [pd_value] = allo_pd_evaluation_final(jacobian_entry, s, s_val, g00, g00_root_val, g01, g01_root_val, g11, g11_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, s, s_val);
    pd_value = subs(pd_value, g00, g00_root_val);
    pd_value = subs(pd_value, g01, g01_root_val);
    pd_value = subs(pd_value, g11, g11_root_val);
end

function [fixed_pt_stabilities] = allo_linear_stability_analysis(jacobian_matrix, s, s_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g00_root_vals;

    jacobian_eval = [0, 0, 0; 
                     0, 0, 0;
                     0, 0, 0];

    for i = 1:length(g00_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = allo_pd_evaluation_final(jacobian_matrix(j, k), s, s_val, g00, g00_root_vals(i), g01, g01_root_vals(i), g11, g11_root_vals(i)); 
            end
        end

        %calulating the eigenvalues of the evaluated jacobian
        eigenvalues = eig(jacobian_eval)

        %classifies the fixed point according to the trace and determinant
        if eigenvalues(1) < 0 && eigenvalues(2) < 0 && eigenvalues(3) < 0
            fixed_pt_stabilities(i) = 1; %1 for stable node
        elseif eigenvalues(1) > 0 && eigenvalues(2) > 0 && eigenvalues(3) > 0
            disp('Error. Unexpected stability type from linear stability analysis.')
        else
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        end
    end
end
