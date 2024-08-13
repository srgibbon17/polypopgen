mu_val = 5e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate

syms s q G0 G1 G2 g0 g1 h mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
%assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1;
w1 = (1-s*h);
w2 = (1-s);

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0*w_bar;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1*w_bar;
%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

mut_g0_coeff = coeffs(simplify(expand(mut_g0)), g0);

mut_g0 = mut_g0/mut_g0_coeff(4);

mut_g0_sub = subs(mut_g0, g0, (g0-(mut_g0_coeff(3)/3)));

mut_g0_sub = (27*s*(2*h - 1))*simplify(expand(mut_g0_sub));

mut_g0_sub_coeff = coeffs(mut_g0_sub, g0);

singularity_set = diff(mut_g0_sub, g0);

isolated_g0 = isolate(singularity_set==0, g0);

equilibria_set = subs(mut_g0_sub, lhs(isolated_g0), rhs(isolated_g0));

% singularity_set = diff(mut_g0, g0);
% 
% isolated_g0 = isolate(singularity_set==0, g0);
% 
% equilibria_set = subs(mut_g0, lhs(isolated_g0), rhs(isolated_g0));

eval_eq_set = subs(equilibria_set, mu, mu_val);
eval_eq_set = subs(eval_eq_set, nu, nu_val);

eval_eq_set = subs(eval_eq_set, h, 1);

s_bifurcations = vpasolve(eval_eq_set==0, s, 0.000001)
s_bifurcations_2 = vpasolve(eval_eq_set==0, s, .99)



