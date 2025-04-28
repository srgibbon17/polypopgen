s_val = 0;
h_val = 0;
mu_val = 0;
nu_val = 0;

syms s q G0 G1 G2 g0 g1 h mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1 & 'Real');
assume(g1>=0 & g1<=1);
assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

%%% g0 is p and g1 is q

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1/w_bar;
w1 = (1-s*h)/w_bar;
w2 = (1-s)/w_bar;

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;

%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

%derivative for linear stability analysis
g0_deriv = diff(mut_g0, g0);


g0_deriv = subs(g0_deriv, mu, mu_val);
g0_deriv = subs(g0_deriv, nu, nu_val);
g0_deriv = subs(g0_deriv, s, s_val);
g0_deriv = subs(g0_deriv, h, h_val);

latex(g0_deriv);

