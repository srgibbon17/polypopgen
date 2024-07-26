% for diploids, a nonlinear model and analysis

% determines number of values to include in 
iterations = 25;


%%% single point parameter values
s_val = 2e-7; % constant value of selection coefficient
mu_val = 1e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate
h_val = 1; % constant value of dominance coefficient
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate

%%% parameter ranges 
s_val_range = logspace(-7, -5, iterations); % range of selection coefficients
mu_val_range = linspace(1e-9, 1e-7, iterations); % range of forward mutation rates
nu_val_range = linspace(1e-9, 1e-7, iterations); % range of backward mutation rates
h_val_range = linspace(.6, 1, iterations); % range of dominance coefficients

syms s q G0 G1 G2 g0 g1 h mu nu

specified_parameter = mu;

variable_list = [g0, specified_parameter];
parameter_list = [s, mu, nu, h];
parameter_values = [s_val, mu_val, nu_val, h_val];

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



G = fsolve_bifn_diagram_init(mut_g0, specified_parameter, parameter_list, parameter_values, g0);
modG = @(x)G(x(1), x(2));

x0_guess = [.5, .000001];
x1_guess = [.01, .000001];
            
options = optimoptions('fsolve', Display='iter');

[x0_soln x0_fval x0_exitflag] = fsolve(modG, x0_guess, options)
[x1_soln x1_fval x1_exitflag] = fsolve(modG, x1_guess, options)

function [diff_eqn_value] = diff_eqn_eval_sym(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0, g0_sub_value)

    diff_eqn_value = subs(mut_exp_g0, mu, mu_value);
    diff_eqn_value = subs(diff_eqn_value, nu, nu_value);
    diff_eqn_value = subs(diff_eqn_value, s, s_value);
    diff_eqn_value = subs(diff_eqn_value, h, h_value);
    diff_eqn_value = subs(diff_eqn_value, g0, g0_sub_value);
end

function [soln] = fixed_points_vpasolve(mut_exp_g0, mu, mu_value, nu, nu_value, s, s_value, h, h_value, g0)

    sub_eqn = subs(mut_exp_g0, mu, mu_value);
    sub_eqn = subs(sub_eqn, nu, nu_value);
    sub_eqn = subs(sub_eqn, s, s_value);
    sub_eqn = subs(sub_eqn, h, h_value);
    
    soln = vpasolve(sub_eqn, g0);
end

function [g0_soln, s_soln] = bifn_vpasolve(delta_g0, derivative_g0, mu, mu_value, nu, nu_value, h, h_value, g0, s)

    delta_g0 = subs(delta_g0, mu, mu_value);
    delta_g0 = subs(delta_g0, nu, nu_value);
    delta_g0 = subs(delta_g0, h, h_value);

    derivative_g0 = subs(derivative_g0, mu, mu_value);
    derivative_g0 = subs(derivative_g0, nu, nu_value);
    derivative_g0 = subs(derivative_g0, h, h_value);
    
    [g0_soln, s_soln] = vpasolve([delta_g0, derivative_g0], [g0, s]);
end

function [delta_g0_value] = delta_g0_eval_implicit(g, s, h, m, n)
    a = s*(1-2*h);
    b = s*(3*h-2-n-h*m+h*n);
    c = -m-n+s-h*s+2*n*s+h*m*s-h*n*s;
    d = n*(1-s);

    delta_g0_value = a*g^3 + b*g^2 + c*g + d;
end



function bifn_surface_plot_pos_sel(G, mu_val_range, h_val_range, nu_val)

    g0_max_array = zeros(iterations, iterations);
    g0_min_array = zeros(iterations, iterations);

    s_max_array = ones(iterations, iterations);
    s_min_array = zeros(iterations, iterations);

    [mu_coord, h_coord] = meshgrid(mu_val_range, h_val_range);

    x0_guess = [.5, .000001];
    x1_guess = [.01, .0000001];

    for i = 1:length(mu_val)
        for j = 1:length(h_val)
            
            [x0_soln x0_fval x0_exitflag] = fsolve(@(x1, x2)G, x0_guess);
            [x1_soln x1_fval x1_exitflag] = fsolve(G, x1_guess);

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

end


function [G] = fsolve_bifn_diagram_init(diff_eqn, specified_parameter, param_list, param_values, g0) 

    syms x1 x2

    %simplifying the difference equation
    diff_eqn = simplify(expand(diff_eqn));

    %calculating symbolic derivative for linear stability analysis
    derivative_sym = diff(diff_eqn, g0);

    %substituting in new, non-symbolic variable into the
    %difference equation and its derivative
    diff_eqn_sub = subs(diff_eqn, g0, x1);
    diff_eqn_sub = subs(diff_eqn_sub, specified_parameter, x2);
    diff_eqn_sub = subs(diff_eqn_sub, param_list(1), param_values(1));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(2), param_values(2));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(3), param_values(3));
    diff_eqn_sub = subs(diff_eqn_sub, param_list(4), param_values(4));

    derivative_sub = subs(derivative_sym, g0, x1);
    derivative_sub = subs(derivative_sub, specified_parameter, x2);
    derivative_sub = subs(derivative_sub, param_list(1), param_values(1));
    derivative_sub = subs(derivative_sub, param_list(2), param_values(2));
    derivative_sub = subs(derivative_sub, param_list(3), param_values(3));
    derivative_sub = subs(derivative_sub, param_list(4), param_values(4));

    F = [1e10*diff_eqn_sub, 1e10*derivative_sub];

    G = matlabFunction(F);
    

end