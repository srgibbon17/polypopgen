% for diploids, a nonlinear model and analysis

iterations = 1000;

s_val = 2e-7; % starting s value

mu_val = 5e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = 1; % h1 dominance coefficient value, constant

syms s q G0 G1 G2 g0 g1 h mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

%%% g0 is q and g1 is p

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1/w_bar;
w1 = (1-s*h)/w_bar;
w2 = (1-s)/w_bar;

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu;
mut_g1 = sel_g0*mu + sel_g1*(1-nu);

% substitution to remove g1
mut_g0 = subs(mut_g0, g1, 1-g0);

g0_values = linspace(0, 1, iterations);
g0_diff_values = [1, iterations];

for i = 1:iterations
    g0_diff_values(i) = diff_eqn_eval(mut_g0, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0, g0_values(i));
end

[g0_soln] = numeric_solver(mut_g0==g0, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0);

h_str = strcat('h: ',string(h_val));
s_str = strcat('s: ',string(s_val));
mu_str = strcat('mu: ',string(mu_val));
nu_str = strcat('nu: ',string(nu_val));
mut_ratio_str = strcat(['mut-ratio ', newline, ...
    '(mu/nu): '],string(mut_ratio_val));

parameters_str = {'Parameters:', s_str, mu_str, nu_str, mut_ratio_str, h_str};
dim = [0.5 0.5 0.3 0.3];

figure
plot(g0_values, g0_diff_values-g0_values, 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'Δg0 function')
hold on
plot(g0_values, zeros(1, length(g0_values)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'g0=Δg0 (nullcline)')

title('Diploids Phase Plane Diagram', 'FontSize', 16)
xlabel('g0', 'FontSize', 14)
ylabel('Δg0/generation', 'FontSize', 14)

plot(g0_soln, zeros(1, length(g0_soln)), 'o', 'Color', "#D95319", 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Fixed Points')


annotation('textbox', dim, 'String', parameters_str,'FontSize', 10, 'FitBoxToText','on')
legend('FontSize', 10)


function [diff_eqn_value] = diff_eqn_eval(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0, g0_sub_value)

    diff_eqn_value = subs(mut_exp_g0, mu, mu_value);
    diff_eqn_value = subs(diff_eqn_value, nu, nu_value);
    diff_eqn_value = subs(diff_eqn_value, s, s_value);
    diff_eqn_value = subs(diff_eqn_value, h, h_value);
    diff_eqn_value = subs(diff_eqn_value, g0, g0_sub_value);
end

function [soln] = numeric_solver(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0)

    sub_eqn = subs(mut_exp_g0, mu, mu_value);
    sub_eqn = subs(sub_eqn, nu, nu_value);
    sub_eqn = subs(sub_eqn, s, s_value);
    sub_eqn = subs(sub_eqn, h, h_value);
    
    soln = vpasolve(sub_eqn, g0);
end
