s_val = 5e-4;
h1_val = 1;
h2_val = 1;
h3_val = 1;
mu_val = 0;
nu_val = 0;
eta_val = 0.001;
kappa_val = 0;


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
eta = sym('eta'); % proportion of HEs in the population
kappa = sym('kappa'); % HE interference/synergy coefficient

% check to ensure eta and kappa are within plausible ranges
if eta_val > 1 || eta_val < 0
    disp('eta_val invalid. eta must be between 0 and 1 (inclusive).');
    return;
elseif kappa_val > 1
    disp('kappa_val invalid. kappa must be less than or equal to 1.');
    return;
elseif kappa_val < -(min(eta_val, 1-eta_val)/max(eta_val, 1-eta_val))
    disp('kappa_val invalid given eta_val. kappa must be greater than or equal to - min(eta, 1-eta)/max(eta, 1-eta).');
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

% equations to parameterize relative fitnesses
wbar = (1-2*s*(h1*(g00*g10+g00*g01)+h2*(g00*g11+g01*g10)+h3*(g01*g11+g10*g11))-s*(h2*(g01^2+g10^2)+g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g00 = g00^2*w0 + (.5+(eta)/16)*2*g01*g00*w1 +   (.5+(eta)/16)*2*g10*g00*w1 +     ((2*eta + eta*(1-eta)*(1-kappa))/8)*g01^2*w2 +     ((2*eta + eta*(1-eta)*(1-kappa))/8)*g10^2*w2 + (.25-(eta*(1-eta)*(1-kappa))/16)*(2*g00*g11+2*g01*g10)*w2 +      ((eta)/16)*2*g01*g11*w3 +      ((eta)/16)*2*g10*g11*w3;
sel_g10 =             ((3*eta)/16)*2*g01*g00*w1 + (.5-(5*eta)/16)*2*g10*g00*w1 + ((eta + (1-kappa)*eta^2 + eta*kappa)/8)*g01^2*w2 +   (1-(6*eta + eta*(1-eta)*(1-kappa))/8)*g10^2*w2 + (.25+(eta*(1-eta)*(1-kappa))/16)*(2*g00*g11+2*g01*g10)*w2 +    ((3*eta)/16)*2*g01*g11*w3 + (.5-(5*eta)/16)*2*g10*g11*w3;
sel_g01 =          (.5-(5*eta)/16)*2*g01*g00*w1 +    ((3*eta)/16)*2*g10*g00*w1 +   (1-(6*eta + eta*(1-eta)*(1-kappa))/8)*g01^2*w2 + ((eta + (1-kappa)*eta^2 + eta*kappa)/8)*g10^2*w2 + (.25+(eta*(1-eta)*(1-kappa))/16)*(2*g00*g11+2*g01*g10)*w2 + (.5-(5*eta)/16)*2*g01*g11*w3 +    ((3*eta)/16)*2*g10*g11*w3;
sel_g11 =               ((eta)/16)*2*g01*g00*w1 +      ((eta)/16)*2*g10*g00*w1 +     ((2*eta + eta*(1-eta)*(1-kappa))/8)*g01^2*w2 +     ((2*eta + eta*(1-eta)*(1-kappa))/8)*g10^2*w2 + (.25-(eta*(1-eta)*(1-kappa))/16)*(2*g00*g11+2*g01*g10)*w2 +   (.5+(eta)/16)*2*g01*g11*w3 +   (.5+(eta)/16)*2*g10*g11*w3 + g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g01*mu*nu + sel_g10*(1-mu)*(1-nu) + sel_g11*(1-nu)*nu - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;

mut_exp_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_exp_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_exp_set(i) = subs(mut_exp_set(i), g10, 1-(g00+g01+g11)); 
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_exp_set(1), g00), diff(mut_exp_set(1), g01), diff(mut_exp_set(1), g11); 
                diff(mut_exp_set(2), g00), diff(mut_exp_set(2), g01), diff(mut_exp_set(2), g11); 
                diff(mut_exp_set(4), g00), diff(mut_exp_set(4), g01), diff(mut_exp_set(4), g11)];


jac_matrix_eval = jac_matrix;

% subs parameter values
for i = 1:length(jac_matrix_eval)
    for j = 1:length(jac_matrix_eval)
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), s, s_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), h1, h1_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), h2, h2_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), h3, h3_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), mu, mu_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), nu, nu_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), eta, eta_val);
        jac_matrix_eval(i,j) = subs(jac_matrix_eval(i,j), kappa, kappa_val);
    end
end

eig_partial_eval = eig(jac_matrix_eval);

% g00 is 10th column followed by g01, g10, g11
sim_data = readmatrix("allo_dominant_fixed.e1e3.k0.txt");

eig1 = zeros(1, 100000);
eig2 = zeros(1, 100000);
eig3 = zeros(1, 100000);

for i = 1:100000
    eig1_val = subs(eig_partial_eval(1), g00, sim_data(i, 10));
    eig1_val = subs(eig1_val, g01, sim_data(i, 11));
    eig1_val = subs(eig1_val, g11, sim_data(i, 13));

    eig1(i) = eig1_val;

    eig2_val = subs(eig_partial_eval(2), g00, sim_data(i, 10));
    eig2_val = subs(eig2_val, g01, sim_data(i, 11));
    eig2_val = subs(eig2_val, g11, sim_data(i, 13));

    eig2(i) = eig2_val;

    eig3_val = subs(eig_partial_eval(3), g00, sim_data(i, 10));
    eig3_val = subs(eig3_val, g01, sim_data(i, 11));
    eig3_val = subs(eig3_val, g11, sim_data(i, 13));

    eig3(i) = eig3_val;

    if mod(i, 1000) == 0
        disp(i);
    end
end

