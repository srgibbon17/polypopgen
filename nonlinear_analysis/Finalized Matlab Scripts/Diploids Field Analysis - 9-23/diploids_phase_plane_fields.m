% for diploids, a nonlinear model and analysis

s_init_val = 7.5e-8;

s_val = s_init_val;
mu_val = 2e-8; % constant value of forward mutation rate
nu_val = 0; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
h_val = 3/4; % h1 dominance coefficient value, constant

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
delta_g1 = subs(mut_g1, g0, 1-g1);

%derivative for linear stability analysis
g1_deriv = diff(delta_g1, g1);

x_index = linspace(0, 1, 100);

figure

%subplot(2, 3, 1)

[g1_roots] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1);

[fixed_pt_stabilities] = linear_stability_analysis(g1_deriv, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1, g1_roots);

plot(x_index, zeros(1, length(x_index)), 'Color','k', 'LineStyle','--', 'LineWidth', 1)
hold on

for i = 1:length(g1_roots)
    if imag(g1_roots(i)) == 0
    if fixed_pt_stabilities(i) == 0
        plot(g1_roots(i), 0, "Color", 'k', 'Marker', 'x', 'LineWidth', 1.25, 'MarkerSize', 12)
    elseif fixed_pt_stabilities(i) == 1
        plot(g1_roots(i), 0, 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 6)
    end
    end
end
    
delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "#0072BD", 'LineWidth', 1)

xlabel 'g_1'
ylabel 'Δ g_1'

s_val = 0;
delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "#D95319", 'LineWidth', 1)

s_val = s_init_val;
mu_val = 0;
nu_val = 0;
delta_g1_values = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, x_index, g1);

plot(x_index, delta_g1_values, 'Color', "#77AC30", 'LineWidth', 1)

xlim([0, 1])





function [delta_g1_values] = delta_g1_eval(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1_index, g1)

    delta_g1 = subs(delta_g1, mu, mu_val);
    delta_g1 = subs(delta_g1, nu, nu_val);
    delta_g1 = subs(delta_g1, s, s_val);
    delta_g1 = subs(delta_g1, h, h_val);

    delta_g1_values = zeros(1, length(g1_index));

    for i = 1:length(g1_index)
        delta_g1_values(i) = subs(delta_g1, g1, g1_index(i));
    end

end

function [g1_root] = root_solns(delta_g1, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g1)

    delta_g1 = subs(delta_g1, mu, mu_val);
    delta_g1 = subs(delta_g1, nu, nu_val);
    delta_g1 = subs(delta_g1, s, s_val);
    delta_g1 = subs(delta_g1, h, h_val);
    
    g1_root = vpasolve(delta_g1, g1);
end

function [fixed_pt_stabilities] = linear_stability_analysis(g0_derivative, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0, g0_roots)
    
    g0_derivative = subs(g0_derivative, mu, mu_val);
    g0_derivative = subs(g0_derivative, nu, nu_val);
    g0_derivative = subs(g0_derivative, s, s_val);
    g0_derivative = subs(g0_derivative, h, h_val);
    
    fixed_pt_stabilities = g0_roots;

    for i = 1:length(fixed_pt_stabilities)
        g0_derivative_eval = subs(g0_derivative, g0, g0_roots(i));
        if g0_derivative_eval < 0
            fixed_pt_stabilities(i) = 1; %1 indicates a stable node
        elseif g0_derivative_eval > 0
            fixed_pt_stabilities(i) = 0; %0 indicates an unstable node
        else
            disp('Error. Linear stabilitiy analysis failed. Derivative is 0.')
        end
    end
end