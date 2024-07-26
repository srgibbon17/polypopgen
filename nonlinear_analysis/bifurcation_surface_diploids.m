% for diploids, a nonlinear model and analysis

iterations = 25;

%s_val = logspace(-7, -5, 1000);
s_val = 2e-7; % starting s value


%mu_val = 1e-8; % constant value of forward mutation rate
mu_val = linspace(1e-9, 1e-7, iterations);
%nu_val = 1e-10; % constant value of backward mutation rate
nu_val = 1e-9;
%mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = linspace(.6, 1, iterations); % h1 dominance coefficient value, constant

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
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;
%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

%derivative for linear stability analysis
derivative_g0 = diff(mut_g0, g0);

g0_max_array = zeros(iterations, iterations);
g0_min_array = zeros(iterations, iterations);

s_max_array = ones(iterations, iterations);
s_min_array = zeros(iterations, iterations);

[mu_coord, h_coord] = meshgrid(mu_val, h_val);

% for i = 1:length(mu_val)
%     for j = 1:length(h_val)
%         [g0_bifn_soln, s_bifn_soln] = bifn_numeric_solver(mut_g0==0, derivative_g0==0, mu, mu_coord(i, j), nu, nu_val, h, h_coord(i, j), g0, s);
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

x0 = [.5, .000001];
x1 = [.01, .0000001];

for i = 1:length(mu_val)
    for j = 1:length(h_val)
        fun_1 = @(x)parameterized_bifn_functions(x, h_coord(i, j), mu_coord(i, j), nu_val);

        [x0_soln x0_fval x0_exitflag] = fsolve(fun_1, x0);
        [x1_soln x1_fval x1_exitflag] = fsolve(fun_1, x1);

        if x0_exitflag > 0
            if x0_soln(1) > 0 && x0_soln(1) < 1 && x0_soln(2) > 0 && x0_soln(2) < 1
                g0_max_array(i, j) = x0_soln(1); 
                s_max_array(i, j) = x0_soln(2);
            end
        end
        if x1_exitflag > 0
            if x1_soln(1) > 0 && x1_soln(1) < 1 && x1_soln(2) > 0 && x1_soln(2) < 1
                g0_min_array(i, j) = x1_soln(1);
                s_min_array(i, j) = x1_soln(2);
            end
        end
    end
end


figure
surf(mu_coord./nu_val, h_coord, s_min_array./s_max_array, abs(g0_min_array-g0_max_array))


xlabel('mu/nu ratio (forward/back)')
ylabel('h (dominance coefficient)')
zlabel('s1/s2 (ratio of selection coefficients)')

title("Diploid Bifurcation Surface")
colorbar


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

function F = parameterized_bifn_functions(x, h, m, n)
    a = x(2)*(1-2*h);
    b = x(2)*(3*h-2-n-h*m+h*n);
    c = -m-n+x(2)-h*x(2)+2*n*x(2)+h*m*x(2)-h*n*x(2);
    d = n*(1-x(2));

    F(1) = 1e10*(a*x(1)^3 + b*x(1)^2 + c*x(1) + d);

    F(2) = 1e10*(3*a*x(1)^2 + 2*b*x(1) + c);

end
