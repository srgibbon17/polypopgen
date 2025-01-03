function [stable_data_matrix] = generalized_allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, beta_val, gamma_val)
%Generates 0HE bifurcation plotting data
%   for allos with 0 HEs, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors

g00 = sym('g00');
g01 = sym('g01');
g10 = sym('g10');
g11 = sym('g11');
s = sym('s');
h1 = sym('h1');
h2 = sym('h2');
h3 = sym('h3');
mu = sym('mu');
nu = sym('nu');
beta = sym('beta');
gamma = sym('gamma');

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
sel_g00 = g00^2*w0 + (.5+(beta+beta*gamma)/32)*2*g01*g00*w1 + (.5+(beta+beta*gamma)/32)*2*g10*g00*w1 + ((3*beta+beta*gamma)/16)*g01^2*w2 + ((3*beta+beta*gamma)/16)*g10^2*w2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + ((beta+beta*gamma)/32)*2*g01*g11*w3 + ((beta+beta*gamma)/32)*2*g10*g11*w3;
sel_g10 = ((3*beta+3*beta*gamma)/32)*2*g01*g00*w1 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g00*w1 + ((beta+3*beta*gamma)/16)*g01^2*w2 + (1-(7*beta+5*beta*gamma)/16)*g10^2*w2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + ((3*beta+3*beta*gamma)/32)*2*g01*g11*w3 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g11*w3;
sel_g01 = (.5-(5*beta+5*beta*gamma)/32)*2*g01*g00*w1 + ((3*beta+3*beta*gamma)/32)*2*g10*g00*w1 + (1-(7*beta+5*beta*gamma)/16)*g01^2*w2 + ((beta+3*beta*gamma)/16)*g10^2*w2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + (.5-(5*beta+5*beta*gamma)/32)*2*g01*g11*w3 + ((3*beta+3*beta*gamma)/32)*2*g10*g11*w3;
sel_g11 = ((beta+beta*gamma)/32)*2*g01*g00*w1 + ((beta+beta*gamma)/32)*2*g10*g00*w1 + ((3*beta+beta*gamma)/16)*g01^2*w2 + ((3*beta+beta*gamma)/16)*g10^2*w2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + (.5+(beta+beta*gamma)/32)*2*g01*g11*w3 + (.5+(beta+beta*gamma)/32)*2*g10*g11*w3 + g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g01*mu*nu + sel_g10*(1-mu)*(1-nu) + sel_g11*(1-nu)*nu - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;

mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g10, 1-(g00+g01+g11));    
end

%creates the Jacobian of the system
jac_matrix = [diff(mut_eqn_set(1), g00), diff(mut_eqn_set(1), g01), diff(mut_eqn_set(1), g11); 
                diff(mut_eqn_set(2), g00), diff(mut_eqn_set(2), g01), diff(mut_eqn_set(2), g11); 
                diff(mut_eqn_set(4), g00), diff(mut_eqn_set(4), g01), diff(mut_eqn_set(4), g11)];

mut_eqn_set_eval = mut_eqn_set;

for i = 1:length(mut_eqn_set_eval)
    % substitutes all constant parameter values (i.e. h1, h2, h3, beta,
    % gamma, mu, and nu)
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), h1, h1_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), h2, h2_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), h3, h3_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), beta, beta_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), gamma, gamma_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), mu, mu_val);
    mut_eqn_set_eval(i) = subs(mut_eqn_set_eval(i), nu, nu_val);
end

jac_matrix_eval = jac_matrix;

for j = 1:length(jac_matrix_eval)
    for k = 1:length(jac_matrix_eval)
        jac_matrix_eval(j,k) = generalized_allo_pd_evaluation_1(jac_matrix_eval(j,k), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val);
    end
end

neutral_g00 = [];
neutral_g01 = [];
neutral_g11 = [];
neutral_s = [];

stable_g00 = [];
stable_g01 = [];
stable_g11 = [];
stable_s = [];

unstable_g00 = [];
unstable_g01 = [];
unstable_g11 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g00_root_vals, g01_root_vals, g11_root_vals] = generalized_allo_root_solns(mut_eqn_set_eval(1), mut_eqn_set_eval(2), mut_eqn_set_eval(4), s, s_val_range(i), g00, g01, g11);
    
    stable_g00(end+1) = g00_root_vals;
    stable_g01(end+1) = g01_root_vals;
    stable_g11(end+1) = g11_root_vals;
    stable_s(end+1) = s_val_range(i);
    

    % %evaluating the jacobian and stability of each fixed point
    % [fixed_pt_stabilities] = generalized_allo_linear_stability_analysis(jac_matrix_eval, s, s_val_range(i), g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);
    % 
    % for j = 1:length(fixed_pt_stabilities)
    % 
    %     if fixed_pt_stabilities(j) == 0
    %         unstable_g00(end+1) = g00_root_vals(j);
    %         unstable_g01(end+1) = g01_root_vals(j);
    %         unstable_g11(end+1) = g11_root_vals(j);
    %         unstable_s(end+1) = s_val_range(i);
    % 
    %     elseif fixed_pt_stabilities(j) == 1
    %         if g00_root_vals(j) > .3333
    %             selected_g00(end+1) = g00_root_vals(j);
    %             selected_g01(end+1) = g01_root_vals(j);
    %             selected_g11(end+1) = g11_root_vals(j);
    %             selected_s(end+1) = s_val_range(i);
    %         else
    %             neutral_g00(end+1) = g00_root_vals(j);
    %             neutral_g01(end+1) = g01_root_vals(j);
    %             neutral_g11(end+1) = g11_root_vals(j);
    %             neutral_s(end+1) = s_val_range(i);
    %         end
    %     end
    % end
end

%neutral_q = neutral_g11 + neutral_g01;
stable_q = stable_g11 + stable_g01;
%unstable_q = unstable_g11 + unstable_g01;

%neutral_g10 = 1 - neutral_g00 - neutral_g01 - neutral_g11;
stable_g10 = 1 - stable_g00 - stable_g01 - stable_g11;
%unstable_g10 = 1 - unstable_g00 - unstable_g01 - unstable_g11;

%neutral_G1 = 2*neutral_g00.*neutral_g10 + 2*neutral_g00.*neutral_g01;
%neutral_G2 = neutral_g10.^2 + 2*(neutral_g00.*neutral_g11 + neutral_g01.*neutral_g10) + neutral_g01.^2;
%neutral_G3 = 2*neutral_g11.*neutral_g10 + 2*neutral_g11.*neutral_g01;
%neutral_genotypes = [neutral_g00.^2; neutral_G1; neutral_G2; neutral_G3; neutral_g11.^2];

stable_G1 = 2*stable_g00.*stable_g10 + 2*stable_g00.*stable_g01;
stable_G2 = stable_g10.^2 + 2*(stable_g00.*stable_g11 + stable_g01.*stable_g10) + stable_g01.^2;
stable_G3 = 2*stable_g11.*stable_g10 + 2*stable_g11.*stable_g01;
stable_genotypes = [stable_g00.^2; stable_G1; stable_G2; stable_G3; stable_g11.^2];

%unstable_G1 = 2*unstable_g00.*unstable_g10 + 2*unstable_g00.*unstable_g01;
%unstable_G2 = unstable_g10.^2 + 2*(unstable_g00.*unstable_g11 + unstable_g01.*unstable_g10) + unstable_g01.^2;
%unstable_G3 = 2*unstable_g11.*unstable_g10 + 2*unstable_g11.*unstable_g01;
%unstable_genotypes = [unstable_g00.^2; unstable_G1; unstable_G2; unstable_G3; unstable_g11.^2];

%neutral_abs_fitness = [ones(1, length(neutral_s)); 1-h1_val*neutral_s; 1-h2_val*neutral_s; 1-h3_val*neutral_s; 1-neutral_s];
stable_abs_fitness = [ones(1, length(stable_s)); 1-h1_val*stable_s; 1-h2_val*stable_s; 1-h3_val*stable_s; 1-stable_s];
%unstable_abs_fitness = [ones(1, length(unstable_s)); 1-h1_val*unstable_s; 1-h2_val*unstable_s; 1-h3_val*unstable_s; 1-unstable_s];

%neutral_avg_fitness = zeros(1, length(neutral_s));
stable_avg_fitness = zeros(1, length(stable_s));
%unstable_avg_fitness = zeros(1, length(unstable_s));

% for i = 1:length(neutral_avg_fitness)
%     neutral_avg_fitness(i) = dot(neutral_genotypes(:, i), neutral_abs_fitness(:, i));
% end

for i = 1:length(stable_avg_fitness)
    stable_avg_fitness(i) = dot(stable_genotypes(:, i), stable_abs_fitness(:, i));
end

% for i = 1:length(unstable_avg_fitness)
%     unstable_avg_fitness(i) = dot(unstable_genotypes(:, i), unstable_abs_fitness(:, i));
% end

stable_data_matrix_transposed = vertcat(stable_s, stable_q, stable_g00, stable_g01, stable_g10, stable_g11, stable_avg_fitness);

stable_data_matrix = stable_data_matrix_transposed.';

end


%%% FUNCTIONS %%%

function [g00_root_vals, g01_root_vals, g11_root_vals] = generalized_allo_root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, s, s_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, s, s_val);
    g01_eqn = subs(mut_g01_eqn, s, s_val);
    g11_eqn = subs(mut_g11_eqn, s, s_val);

    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end


function [pd_value] = generalized_allo_pd_evaluation_1(jacobian_entry, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val)
    
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

function [pd_value] = generalized_allo_pd_evaluation_2(jacobian_entry, s, s_val, g00, g00_root_val, g01, g01_root_val, g11, g11_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, s, s_val);
    pd_value = subs(pd_value, g00, g00_root_val);
    pd_value = subs(pd_value, g01, g01_root_val);
    pd_value = subs(pd_value, g11, g11_root_val);
end

function [fixed_pt_stabilities] = generalized_allo_linear_stability_analysis(jacobian_matrix, s, s_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals)
    
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
                jacobian_eval(j, k) = generalized_allo_pd_evaluation_2(jacobian_matrix(j, k), s, s_val, g00, g00_root_vals(i), g01, g01_root_vals(i), g11, g11_root_vals(i)); 
            end
        end

        %calulating the eigenvalues of the evaluated jacobian
        [eigenvectors, eigenvalues] = eig(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if eigenvalues(1,1) < 0 && eigenvalues(2,2) < 0 && eigenvalues(3,3) < 0
            fixed_pt_stabilities(i) = 1; %0 for stable node
        elseif eigenvalues(1,1) > 0 && eigenvalues(2,2) > 0 && eigenvalues(3,3) > 0
            disp('Error. Unexpected stability type from linear stability analysis.')
        else
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        end
    end

end