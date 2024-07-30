% for diploids, a nonlinear model and analysis

iterations = 100;

%s_val = linspace(1e-8, 1e-5, 100);
s_val = 1e-7; % starting s value

mu_val = 1e-8; % constant value of forward mutation rate
%mu_val = logspace(-9, -7, 100);
nu_val = 1e-10; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = 1; % h1 dominance coefficient value, constant

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
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;
%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

%derivative for linear stability analysis
derivative_g0 = diff(mut_g0, g0);

%creating range of g0_values to evalute at
g0_values = linspace(0, 1, iterations);
g0_diff_values = [1, iterations];


% for i = 1:iterations
%     g0_diff_values(i) = diff_eqn_eval(mut_g0, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0, g0_values(i));
% end
% 
% figure
% plot(g0_values, g0_diff_values)
% hold on
% plot(g0_values, zeros(1, length(g0_values)), 'LineStyle','--', 'Color', 'k')
% title('g0 Nonlinear Analysis')
% xlabel('g0')
% ylabel('delta g0 over generations')



[g0_bifn_soln, s_bifn_soln] = bifn_numeric_solver(mut_g0==0, derivative_g0==0, nu, nu_val, mu, mu_val, s, s_val, g0, h)

g0_bifn_values = [];
s_bifn_values = [];
for i = 1:length(g0_bifn_soln)
    if g0_bifn_soln(i) ~= 0
        g0_bifn_values(end+1) = g0_bifn_soln(i);
        s_bifn_values(end+1) = s_bifn_soln(i);
    end
end

g0_stable_1 = [];
s_stable_1 = [];

g0_stable_2 = [];
s_stable_2 = [];

g0_unstable = [];
s_unstable = [];

[g0_soln] = numeric_solver(mut_g0==0, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0);

% for i = 1:length(mu_val)
%     for j = 1:length(nu_val)
%         for k = 1:length(s_val)
%             for l = 1:length(h_val)
%                 g0_soln_real = [];
%                 s_soln_real = [];
%                 [g0_soln] = numeric_solver(mut_g0==0, mu, mu_val(i), nu, nu_val(j), s, s_val(k), h, h_val(l), g0);
%                 for m = 1:length(g0_soln)
%                     if imag(g0_soln(m)) == 0
%                         g0_soln_real(end+1) = g0_soln(m);
%                         s_soln_real(end+1) = s_val(k);
%                     end
%                 end
%                 for m = 1:length(g0_soln_real)
%                     if g0_soln_real(m) < min(g0_bifn_values)
%                         g0_stable_1(end+1) = g0_soln_real(m);
%                         s_stable_1(end+1) = s_soln_real(m);
%                     elseif g0_soln_real(m) > min(g0_bifn_values) && g0_soln_real(m) < max(g0_bifn_values)
%                         g0_unstable(end+1) = g0_soln_real(m);
%                         s_unstable(end+1) = s_soln_real(m);
%                     elseif g0_soln_real(m) > max(g0_bifn_values)
%                         g0_stable_2(end+1) = g0_soln_real(m);
%                         s_stable_2(end+1) = s_soln_real(m);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% h_str = strcat('h: ',string(h_val));
% mu_str = strcat('mu: ',string(mu_val));
% nu_str = strcat('nu: ',string(nu_val));
% mut_ratio_str = strcat(['mut-ratio ', newline, ...
%     '(mu/nu): '],string(mut_ratio_val));
% 
% parameters_str = {'Parameters:', mu_str, nu_str, mut_ratio_str, h_str};
% dim = [0.5 0.5 0.3 0.3];
% 
% figure 
% 
% subplot(1, 2, 1)
% 
% plot(s_stable_1, g0_stable_1, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable')
% hold on
% plot(s_unstable, g0_unstable, 'LineStyle','--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Unstable')
% plot(s_stable_2, g0_stable_2, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable')
% plot(s_bifn_values(1), g0_bifn_values(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point 1')
% plot(s_bifn_values(2), g0_bifn_values(2), '.', 'Color', '#0072BD','MarkerSize', 20, 'DisplayName', 'Bifurcation Point 2')
% 
% xlabel('s', 'FontSize', 14)
% ylabel('g0', 'FontSize', 14)
% sgtitle('Diploid Bifurcation Diagram', 'FontSize', 16)  
% xscale log
% 
% %xlim([0 1e-6])
% %ylim([0 1])
% annotation('textbox', dim, 'String', parameters_str,'FontSize', 10, 'FitBoxToText','on')
% 
% subplot(1, 2, 2)
% 
% g1_stable_1 = ones(1, length(g0_stable_1)) - g0_stable_1;
% g1_stable_2 = ones(1, length(g0_stable_2)) - g0_stable_2;
% g1_unstable = ones(1, length(g0_unstable)) - g0_unstable;
% 
% plot(s_stable_1, g1_stable_1, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable')
% hold on
% plot(s_unstable, g1_unstable, 'LineStyle','--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Unstable')
% plot(s_stable_2, g1_stable_2, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'Stable')
% plot(s_bifn_values(1), 1-g0_bifn_values(1), '.', 'Color', '#0072BD', 'MarkerSize', 20, 'DisplayName', 'Bifurcation Point 1')
% plot(s_bifn_values(2), 1-g0_bifn_values(2), '.', 'Color', '#0072BD','MarkerSize', 20, 'DisplayName', 'Bifurcation Point 2')
% 
% xlabel('s', 'FontSize', 14)
% xscale log
% 
% %xlim([0 3e-8])
% %ylim([0 1])


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

function [g0_soln, s_soln] = bifn_numeric_solver(delta_g0, derivative_g0, mu, mu_value, nu, nu_value, h, h_value, g0, s)

    delta_g0 = subs(delta_g0, mu, mu_value);
    delta_g0 = subs(delta_g0, nu, nu_value);
    delta_g0 = subs(delta_g0, h, h_value);

    derivative_g0 = subs(derivative_g0, mu, mu_value);
    derivative_g0 = subs(derivative_g0, nu, nu_value);
    derivative_g0 = subs(derivative_g0, h, h_value);
    
    [g0_soln, s_soln] = vpasolve([delta_g0, derivative_g0], [g0, s]);
end
