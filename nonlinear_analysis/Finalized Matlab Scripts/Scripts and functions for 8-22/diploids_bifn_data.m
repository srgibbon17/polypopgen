function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val)
%Generates diploid bifurcation plotting data
%   for diploids, a nonlinear model and analysis
%   returns:
%       diploid_data:
%           a 6 element array of data arrays in the following order:
%               neutral_stable (for q and s, respectively)
%               selected_stable (for q and s, respectively)
%               unstable (for q and s, respectively)

g0 = sym('g0');
g1 = sym('g1');
G0 = sym('G0');
G1 = sym('G1');
G2 = sym('G2');
s = sym('s');
h = sym('h');
mu = sym('mu');
nu = sym('nu');

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

neutral_stable_g0 = [];
neutral_stable_s = [];

selected_stable_g0 = [];
selected_stable_s = [];

unstable_g0 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    [g0_roots] = dip_root_solns(mut_g0, mu, mu_val, nu, nu_val, s, s_val_range(i), h, h_val, g0);

    [root_stabilities] = dip_linear_stability_analysis(g0_deriv, mu, mu_val, nu, nu_val, s, s_val_range(i), h, h_val, g0, g0_roots);

    for j = 1:length(root_stabilities)
        if root_stabilities(j) == 0
            if imag(g0_roots(j)) == 0 && real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                unstable_g0(end+1) = g0_roots(j);
                unstable_s(end+1) = s_val_range(i);
            end

        elseif root_stabilities(j) == 1
            if imag(g0_roots(j)) == 0 
                if s_val_range(i) > 0 && g0_roots(j) > 1/3 && real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                    selected_stable_g0(end+1) = g0_roots(j);
                    selected_stable_s(end+1) = s_val_range(i);
                elseif s_val_range(i) < 0 && g0_roots(j) > 2/3 && real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                    selected_stable_g0(end+1) = g0_roots(j);
                    selected_stable_s(end+1) = s_val_range(i);
                elseif real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                    neutral_stable_g0(end+1) = g0_roots(j);
                    neutral_stable_s(end+1) = s_val_range(i);
                end
            end
        end
    end
end

neutral_stable_q = ones(1, length(neutral_stable_g0)) - neutral_stable_g0;
selected_stable_q = ones(1, length(selected_stable_g0)) - selected_stable_g0;
unstable_q = ones(1, length(unstable_g0)) - unstable_g0;

end

function [g0_root] = dip_root_solns(mut_eqn_g0, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0)

    mut_eqn_g0 = subs(mut_eqn_g0, mu, mu_val);
    mut_eqn_g0 = subs(mut_eqn_g0, nu, nu_val);
    mut_eqn_g0 = subs(mut_eqn_g0, s, s_val);
    mut_eqn_g0 = subs(mut_eqn_g0, h, h_val);
    
    g0_root = vpasolve(mut_eqn_g0, g0);
end

function [fixed_pt_stabilities] = dip_linear_stability_analysis(g0_derivative, mu, mu_val, nu, nu_val, s, s_val, h, h_val, g0, g0_roots)
    
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