% for autos, classification of fixed points using linear stability
% analysis, the Jacobian matrix, and eigenvectors

s_val = 0;
mu_val = 0; % constant value of forward mutation rate
nu_val = 0; % constant value of backward mutation rate

h1_val = 0; % h1 dominance coefficient value, constant
h2_val = 0; % h2 dominance coefficient value, constant
h3_val = 0; % h3 dominance coefficient value, constant

syms s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu

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
sel_g0 = G0*w0+(1/2)*G1*w1 + (1/6)*G2*w2;
sel_g1 = (1/2)*G1*w1 + (2/3)*G2*w2 + (1/2)*G3*w3;
sel_g2 = (1/6)*G2*w2 + (1/2)*G3*w3 + G4*w4;

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

jac_matrix_eval = jac_matrix;

for i = 1:length(jac_matrix_eval)
    for j = 1:length(jac_matrix_eval)
        jac_matrix_eval(i,j) = pd_evaluation(jac_matrix(i,j), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val);
    end
end

[eigenspaces, eigenvalues] = eig(jac_matrix_eval)

function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val)
    
    %%%function which evaluates a partial derivative by substituting for
    %%%parameters
    %is used to evaluate the jacobian

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, s, s_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
end

