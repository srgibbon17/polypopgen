mu_val = 5e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate

syms s q G0 G1 G2 g0 g1 h mu nu g0p

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
mut_g0 = simplify(expand(subs(mut_g0, g1, 1-g0)));

first_der = diff(mut_g0, g0);
second_der = diff(first_der, g0);

eqn_set = [mut_g0, first_der, second_der];

for i = 1:length(eqn_set)
   eqn_set(i) = subs(eqn_set(i), mu, mu_val); 
   eqn_set(i) = subs(eqn_set(i), nu, nu_val);
end

S = vpasolve(eqn_set, [g0, s, h]);

% coeffcients of g0 from lowest to highest degree (i.e. constant, linear,
% quadratic, cubic)
% mut_g0_coeffs = coeffs(mut_g0, g0);
% 
% % from the cubic formula/algebra, a double root occurs 
% % if 9ad-bc = 2(b^2 - 3ac)
% % for any cubic of the general form a*x^3 + b*x^2 + c*x + d
% % this corresponds to the r paramter > 0 in Strogatz 3.6
% 
% a = mut_g0_coeffs(4);
% b = mut_g0_coeffs(3);
% c = mut_g0_coeffs(2);
% d = mut_g0_coeffs(1);
% 
% bifn_curves = 2*(b^2 - 3*a*c) == 9*a*d - b*c;
% 
% bifn_curves = subs(bifn_curves, mu, mu_val);
% bifn_curves = subs(bifn_curves, nu, nu_val);
% 
% s_bifn_curves = solve(bifn_curves, s);
% 
% discriminant = 18*a*b*c*d - 4*b^3*d + b^2*c^2 - 4*a*c^3 - 27*a^2*d^2;
% 
% bifn_curves_2 = (9*a*d - b*c)/2*(b^2 - 3*a*c);
% 
% bifn_curves_solns = bifn_curves_2 == discriminant;
% 
% bifn_curves_solns = subs(bifn_curves_solns, mu, mu_val);
% bifn_curves_solns = subs(bifn_curves_solns, nu, nu_val);
% 
% triple_root = discriminant == b^2 - 3*a*c;
% 
% triple_root_solns = subs(triple_root, mu, mu_val);
% triple_root_solns = subs(triple_root_solns, nu, nu_val);