syms alpha s q G_0 G_1 G_2 G_3 G_4 g_0 g_1 g_2 h_1 h_2 h_3 mu nu omega

% assumptions on the parameters of the model; theoretical bounds
assume(g_0>=0 & g_0<=1);
assume(g_1>=0 & g_1<=1);
assume(g_2>=0 & g_2<=1);
assume(s>=-1 & s<=1);
assume(h_1>=0 & h_1<=1);
assume(h_2>=0 & h_2<=1);
assume(h_3>=0 & h_3<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(alpha>=0 & alpha<=1/6);
assume(G_0>=0 & G_0<=1);
assume(G_1>=0 & G_1<=1);
assume(G_2>=0 & G_2<=1);
assume(G_3>=0 & G_3<=1);
assume(G_4>=0 & G_4<=1);

% equations to parameterize relative fitnessesbar
%omega = 1 - s*(G_1*h_1 + G_2*h_2 + G_3*h_3 + G_4);
w0 = 1;
w1 = (1-s*h_1);
w2 = (1-s*h_2);
w3 = (1-s*h_3);
w4 = (1-s);

% equations for selection
sel_g0 = G_0*w0+(1/2 + alpha/4)*G_1*w1 + (1/6 + alpha/3)*G_2*w2 + (alpha/4)*G_3*w3;
sel_g1 = (1/2 - alpha/2)*G_1*w1 + (2/3 - 2*alpha/3)*G_2*w2 + (1/2 - alpha/2)*G_3*w3;
sel_g2 = (alpha/4)*G_1*w1 + (1/6 + alpha/3)*G_2*w2 + (1/2 + alpha/4)*G_3*w3 + G_4*w4;

% equations for mutation
mut_g0 = (sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2))*(1/omega) - g_0;
mut_g1 = (2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu)*(1/omega) - g_1;
mut_g2 = (sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2))*(1/omega) - g_2;

mut_exp_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), G_0, g_0^2);
    mut_exp_set(i) = subs(mut_exp_set(i), G_1, 2*g_0*g_1);
    mut_exp_set(i) = subs(mut_exp_set(i), G_2, (2*g_0*g_2 + g_1^2));
    mut_exp_set(i) = subs(mut_exp_set(i), G_3, 2*g_1*g_2);
    mut_exp_set(i) = subs(mut_exp_set(i), G_4, g_2^2);
    mut_exp_set(i) = subs(mut_exp_set(i), g_2, (1-g_1-g_0));
end


%creates the Jacobian of the system
jacobian_matrix = [diff(mut_exp_set(1), g_0), diff(mut_exp_set(1), g_1); 
                   diff(mut_exp_set(2), g_0), diff(mut_exp_set(2), g_1)];

disp('Unsimplified Jacobian First Entry')
disp(jacobian_matrix(1,1))

disp(latex(jacobian_matrix(1,1)))

jacobian_matrix_simp = jacobian_matrix;

for i = 1:2
    for j = 1:2
        % "removes" nu from the expression by setting it equal to zero
        % this line can be copied and nu and can be replaced to remove
        % alpha or dominance coefficients as well
        jacobian_matrix_simp(i, j) = subs(jacobian_matrix_simp(i, j), nu, 0);
    end
end

disp('Simplified Jacobian First Entry (nu and other parameters set to 0)')

disp(jacobian_matrix_simp(1,1))

disp(latex(jacobian_matrix_simp(1,1)))