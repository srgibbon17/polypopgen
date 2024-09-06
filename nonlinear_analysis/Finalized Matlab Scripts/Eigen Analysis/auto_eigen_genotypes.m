function auto_eigen_genotypes(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val)
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
g0_mut = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2);
g1_mut = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu;
g2_mut = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2);

G0_t_1 = (g0_mut)^2;
G1_t_1 = 2*g0_mut*g1_mut;
G2_t_1 = 2*g0_mut*g2_mut + (g1_mut)^2;
G3_t_1 = 2*g1_mut*g2_mut;
G4_t_1 = (g2_mut)^2;

delta_G0 = G0_t_1 - G0;
delta_G1 = G1_t_1 - G1;
delta_G2 = G2_t_1 - G2;
delta_G3 = G3_t_1 - G3;
delta_G4 = G4_t_1 - G4;

diff_eqn_set = [delta_G0, delta_G1, delta_G2, delta_G3, delta_G4];

for i = 1:length(diff_eqn_set)
    diff_eqn_set(i) = subs(diff_eqn_set(i), G2, (1-G0-G1-G3-G4));
end

%jac_matrix = [diff(diff_eqn_set(1), G0), diff(diff_eqn_set(1), G1), diff(diff_eqn_set(1), G3), diff(diff_eqn_set(1), G4);
%              diff(diff_eqn_set(2), G0), diff(diff_eqn_set(2), G1), diff(diff_eqn_set(2), G3), diff(diff_eqn_set(2), G4);
%              diff(diff_eqn_set(4), G0), diff(diff_eqn_set(4), G1), diff(diff_eqn_set(4), G3), diff(diff_eqn_set(4), G4);
%              diff(diff_eqn_set(5), G0), diff(diff_eqn_set(5), G1), diff(diff_eqn_set(5), G3), diff(diff_eqn_set(5), G4)];

% equations for mutation
g0_t_1 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2);
g1_t_1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu;
g2_t_1 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2);

delta_g0 = g0_t_1 - g0;
delta_g1 = g1_t_1 - g1;
delta_g2 = g2_t_1 - g2;

diff_eqn_set_2 = [delta_g0, delta_g1, delta_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(diff_eqn_set_2)
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), G0, g0^2);
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), G1, 2*g0*g1);
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), G2, (2*g0*g2 + g1^2));
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), G3, 2*g1*g2);
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), G4, g2^2);
    diff_eqn_set_2(i) = subs(diff_eqn_set_2(i), g2, (1-g1-g0));
end




for i = 1:length(s_val_range)
    disp(i)
    [g0_root_vals, g1_root_vals] = root_solns_2(diff_eqn_set_2(1), diff_eqn_set_2(2), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1)

    g2_root_vals = 1- g0_root_vals - g1_root_vals;

    initial_guess = [g0_root_vals^2, 2*g0_root_vals*g1_root_vals, 2*g1_root_vals*g2_root_vals, g2_root_vals^2]

    [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = root_solns(diff_eqn_set(1), diff_eqn_set(2), diff_eqn_set(4), diff_eqn_set(5), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, G0, G1, G3, G4, initial_guess)
end

end

function [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = root_solns(diff_G0, diff_G1, diff_G3, diff_G4, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, G0, G1, G3, G4, initial_guess)

    %function which uses vpasolve to find the fixed points/roots of the system

    G0_eqn = subs(diff_G0, mu, mu_val);
    G0_eqn = subs(G0_eqn, nu, nu_val);
    G0_eqn = subs(G0_eqn, s, s_val);
    G0_eqn = subs(G0_eqn, h1, h1_val);
    G0_eqn = subs(G0_eqn, h2, h2_val);
    G0_eqn = subs(G0_eqn, h3, h3_val);
    G0_eqn = subs(G0_eqn, a, a_val);

    G1_eqn = subs(diff_G1, mu, mu_val);
    G1_eqn = subs(G1_eqn, nu, nu_val);
    G1_eqn = subs(G1_eqn, s, s_val);
    G1_eqn = subs(G1_eqn, h1, h1_val);
    G1_eqn = subs(G1_eqn, h2, h2_val);
    G1_eqn = subs(G1_eqn, h3, h3_val);
    G1_eqn = subs(G1_eqn, a, a_val);

    G3_eqn = subs(diff_G3, mu, mu_val);
    G3_eqn = subs(G3_eqn, nu, nu_val);
    G3_eqn = subs(G3_eqn, s, s_val);
    G3_eqn = subs(G3_eqn, h1, h1_val);
    G3_eqn = subs(G3_eqn, h2, h2_val);
    G3_eqn = subs(G3_eqn, h3, h3_val);
    G3_eqn = subs(G3_eqn, a, a_val);

    G4_eqn = subs(diff_G4, mu, mu_val);
    G4_eqn = subs(G4_eqn, nu, nu_val);
    G4_eqn = subs(G4_eqn, s, s_val);
    G4_eqn = subs(G4_eqn, h1, h1_val);
    G4_eqn = subs(G4_eqn, h2, h2_val);
    G4_eqn = subs(G4_eqn, h3, h3_val);
    G4_eqn = subs(G4_eqn, a, a_val);


    [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = vpasolve([G0_eqn, G1_eqn, G3_eqn, G4_eqn], [G0, G1, G3, G4], initial_guess);
end

function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, G0, G0_root_val, G1, G1_root_val, G3, G3_root_val, G4, G4_root_val)
    
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
    pd_value = subs(pd_value, G0, G0_root_val);
    pd_value = subs(pd_value, G1, G1_root_val);
    pd_value = subs(pd_value, G3, G3_root_val);
    pd_value = subs(pd_value, G4, G4_root_val);
end

function [g0_root_vals, g1_root_vals] = root_solns_2(mut_g0_eqn, mut_g1_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1)

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


