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
    %mut_exp_set(i) = subs(mut_exp_set(i), g2, (1-g1-g0));
end

for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), g0, (1-g1-g2));
    mut_exp_set(i) = subs(mut_exp_set(i), g1, 2*(q-g2));
    mut_exp_set(i) = subs(mut_exp_set(i), g2, (q^2));
end

s_value = 1e-5;
mu_value = 5e-9;
nu_value = 1e-10;
h1_value = .05;
h2_value = .7;
h3_value = .9;
a_value = 1/12;

for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), s, s_value);
    mut_exp_set(i) = subs(mut_exp_set(i), mu, mu_value);
    mut_exp_set(i) = subs(mut_exp_set(i), nu, nu_value);
    mut_exp_set(i) = subs(mut_exp_set(i), h1, h1_value);
    mut_exp_set(i) = subs(mut_exp_set(i), h2, h2_value);
    mut_exp_set(i) = subs(mut_exp_set(i), h3, h3_value);
    %mut_exp_set(i) = subs(mut_exp_set(i), a, a_value);
    
end

[q_soln, a_soln] = vpasolve([mut_exp_set(1), mut_exp_set(2), mut_exp_set(3)], [q, a]);

%q_soln = solve(mut_exp_set(1), q, 'ReturnConditions', true);