% for diploids, a nonlinear model and analysis

iterations = 25;

%s_val = logspace(-7, -5, 1000);
s_val = 2e-7; % starting s value

mu_val = linspace(8e-9, 8.002e-9, 50); % constant value of forward mutation rate
%mu_val = logspace(-9, -7, iterations);
nu_val = 1e-9; % constant value of backward mutation rate
%nu_val = logspace(-10, -8, iterations);
%mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = 1; % h1 dominance coefficient value, constant

% syms s q G0 G1 G2 g0 g1 h mu nu
% 
% % assumptions on the parameters of the model; theoretical bounds
% assume(g0>=0 & g0<=1 & 'Real');
% assume(g1>=0 & g1<=1);
% assume(s>=-1 & s<=1);
% assume(h>=0 & h<=1);
% assume(mu>=0 & mu<=1);
% assume(nu>=0 & nu<=1);
% assume(G0>=0 & G0<=1);
% assume(G1>=0 & G1<=1);
% assume(G2>=0 & G2<=1);
% 
% %%% g0 is q and g1 is p
% 
% % equations to parameterize relative fitnesses
% w_bar = 1 - s*(h*2*g0*g1 + g1^2);
% w0 = 1/w_bar;
% w1 = (1-s*h)/w_bar;
% w2 = (1-s)/w_bar;
% 
% % equations for selection
% sel_g0 = w1*g1*g0 + w0*g0^2; 
% sel_g1 = w2*g1^2 + w1*g1*g0;
% 
% % equations for mutation
% mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
% mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;
% %removing g1 from the equations
% mut_g0 = subs(mut_g0, g1, 1-g0);
% 
% %derivative for linear stability analysis
% derivative_g0 = diff(mut_g0, g0);
% 
% for i = 1:length(mu_val)
%     disp(strcat('mu value: ', string(mu_val(i))))
%     [g0_bifn_soln, s_bifn_soln] = bifn_numeric_solver(mut_g0==0, derivative_g0==0, mu, mu_val(i), nu, nu_val, h, h_val, g0, s)
% end

% g0_max_array = ones(iterations, iterations);
% g0_min_array = ones(iterations, iterations);
% 
% s_max_array = ones(iterations, iterations);
% s_min_array = zeros(iterations, iterations);
% 
% [mu_coord, nu_coord] = meshgrid(mu_val, nu_val);
% 
% for i = 1:length(mu_val)
%     for j = 1:length(nu_val)
%         [g0_bifn_soln, s_bifn_soln] = bifn_numeric_solver(mut_g0==0, derivative_g0==0, mu, mu_coord(i, j), nu, nu_coord(i, j), h, h_val, g0, s);
%         g0_bifn_values = [];
%         s_bifn_values = [];
%         if length(g0_bifn_soln) > 1
%             for k = 1:length(g0_bifn_soln)
%                 if g0_bifn_soln(k) ~= 0
%                     g0_bifn_values(end+1) = g0_bifn_soln(k);
%                     s_bifn_values(end+1) = s_bifn_soln(k);    
%                 end
%             end
%             g0_max_array(i, j) = max(g0_bifn_values); 
%             g0_min_array(i, j) = min(g0_bifn_values);
%             for k = 1:length(g0_bifn_values)
%                 if g0_bifn_values(k) == max(g0_bifn_values)
%                     s_max_array(i, j) = s_bifn_values(k);
%                 elseif g0_bifn_values(k) == min(g0_bifn_values)
%                     s_min_array(i,j) = s_bifn_values(k);
%                 else
%                     disp('Error. Unexpected number of bifurcation points.')
%                 end
%             end
%         end
%     end
% end


% 
% for i = 1:length(h_val)
%     disp(h_val(i))
%     [g0_bifn_soln, s_bifn_soln] = bifn_numeric_solver(mut_g0==0, derivative_g0==0, mu, mu_val, nu, nu_val, h, h_val(i), g0, s)
% end


% g0_stable_1 = [];
% s_stable_1 = [];
% 
% g0_stable_2 = [];
% s_stable_2 = [];
% 
% g0_unstable = [];
% s_unstable = [];

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
% xlim([0 1e-6])
% ylim([0 1])
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
% xlim([0 1e-6])
% ylim([0 1])

% figure
% surf(mu_coord, nu_coord, s_min_array./s_max_array, abs(g0_min_array-g0_max_array))
% 
% xlabel('mu (forward mutation)')
% ylabel('nu (back mutation)')
% zlabel('s2/s1')
% xscale log
% yscale log
% 
% title("Diploid Bifurcation Surface")
% colorbar

h = 0;
mu = 1e-9;
nu = 1e-7;
s = 1e-7;
k = 1e10;

fun_1 = @(x)parameterized_bifn_functions(x, h, mu, nu, k);
x0 = [.5, -.000001];
x1 = [.999, -.000001];

options = optimoptions('fsolve', 'Display', 'iter', 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6);

[x0_soln fval exitflag] = fsolve(fun_1, x0)
[x1_soln fval exitflag] = fsolve(fun_1, x1, options)

% fun_2 = @(x)parameterized_functions(x, s, h, mu, nu, k);
% 
% x0_prime = 0.0001;
% x1_prime = .45;
% x2_prime = .999;
% 
% x0_soln = fsolve(fun_2, x0_prime, options)
% x1_soln = fsolve(fun_2, x1_prime, options)
% x2_soln = fsolve(fun_2, x2_prime, options)


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

%%% new functions for fsolve

function F = parameterized_bifn_functions(x, h, m, n, k)
    a = x(2)*(1-2*h);
    b = x(2)*(3*h-2-n-h*m+h*n);
    c = -m-n+x(2)-h*x(2)+2*n*x(2)+h*m*x(2)-h*n*x(2);
    d = n*(1-x(2));

    F(1) = k*(a*x(1)^3 + b*x(1)^2 + c*x(1) + d);

    F(2) = k*(3*a*x(1)^2 + 2*b*x(1) + c);

end

function F = parameterized_functions(x, s, h, m, n, k)
    a = s*(1-2*h);
    b = s*(3*h-2-n-h*m+h*n);
    c = -m-n+s-h*s+2*n*s+h*m*s-h*n*s;
    d = n*(1-s);
    

    F = k*(a*x^3 + b*x^2 + c*x + d);

end
