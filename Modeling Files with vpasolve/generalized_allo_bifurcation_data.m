%Generates bifurcation plotting data for a generalized allo model under
%the beta, gamma parameterization
%   classifies fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
   

g00 = sym('g00');
g01 = sym('g01');
g10 = sym('g10');
g11 = sym('g11');
s = sym('s');
h1 = sym('h1');
h2 = sym('h2');
h3 = sym('h3');
mu = sym('mu');
nu = sym('nu');
beta = sym('beta');
gamma = sym('gamma');

% assumptions on the parameters of the model; theoretical bounds
assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);

% equations to parameterize relative fitnesses
wbar = (1-2*s*(h1*(g00*g10+g00*g01)+h2*(g00*g11+g01*g10)+h3*(g01*g11+g10*g11))-s*(h2*(g01^2+g10^2)+g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g00 = g00^2*w0 + (.5+(beta+beta*gamma)/32)*2*g01*g00*w1 + (.5+(beta+beta*gamma)/32)*2*g10*g00*w1 + ((3*beta+beta*gamma)/16)*g01^2*w2 + ((3*beta+beta*gamma)/16)*g10^2*w2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + ((beta+beta*gamma)/32)*2*g01*g11*w3 + ((beta+beta*gamma)/32)*2*g10*g11*w3;
sel_g10 = ((3*beta+3*beta*gamma)/32)*2*g01*g00*w1 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g00*w1 + ((beta+3*beta*gamma)/16)*g01^2*w2 + (1-(7*beta+5*beta*gamma)/16)*g10^2*w2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + ((3*beta+3*beta*gamma)/32)*2*g01*g11*w3 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g11*w3;
sel_g01 = (.5-(5*beta+5*beta*gamma)/32)*2*g01*g00*w1 + ((3*beta+3*beta*gamma)/32)*2*g10*g00*w1 + (1-(7*beta+5*beta*gamma)/16)*g01^2*w2 + ((beta+3*beta*gamma)/16)*g10^2*w2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + (.5-(5*beta+5*beta*gamma)/32)*2*g01*g11*w3 + ((3*beta+3*beta*gamma)/32)*2*g10*g11*w3;
sel_g11 = ((beta+beta*gamma)/32)*2*g01*g00*w1 + ((beta+beta*gamma)/32)*2*g10*g00*w1 + ((3*beta+beta*gamma)/16)*g01^2*w2 + ((3*beta+beta*gamma)/16)*g10^2*w2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10)*w2 + (.5+(beta+beta*gamma)/32)*2*g01*g11*w3 + (.5+(beta+beta*gamma)/32)*2*g10*g11*w3 + g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g01*mu*nu + sel_g10*(1-mu)*(1-nu) + sel_g11*(1-nu)*nu - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;

mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g10, 1-(g00+g01+g11));    
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_eqn_set(1), g00), diff(mut_eqn_set(1), g01), diff(mut_eqn_set(1), g11); 
                diff(mut_eqn_set(2), g00), diff(mut_eqn_set(2), g01), diff(mut_eqn_set(2), g11); 
                diff(mut_eqn_set(4), g00), diff(mut_eqn_set(4), g01), diff(mut_eqn_set(4), g11)];

neutral_stable_g00 = [];
neutral_stable_g01 = [];
neutral_stable_g11 = [];
neutral_stable_s = [];

selected_stable_g00 = [];
selected_stable_g01 = [];
selected_stable_g11 = [];
selected_stable_s = [];

unstable_g00 = [];
unstable_g01 = [];
unstable_g11 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g00_root_vals, g01_root_vals, g11_root_vals] = generalized_allo_root_solns(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = generalized_allo_linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);

    for j = 1:length(fixed_pt_stabilities)

        if fixed_pt_stabilities(j) == 0
            unstable_g00(end+1) = g00_root_vals(j);
            unstable_g01(end+1) = g01_root_vals(j);
            unstable_g11(end+1) = g11_root_vals(j);
            unstable_s(end+1) = s_val_range(i);

        elseif fixed_pt_stabilities(j) == 1
            if g00_root_vals(j) > .3333
                selected_stable_g00(end+1) = g00_root_vals(j);
                selected_stable_g01(end+1) = g01_root_vals(j);
                selected_stable_g11(end+1) = g11_root_vals(j);
                selected_stable_s(end+1) = s_val_range(i);
            else
                neutral_stable_g00(end+1) = g00_root_vals(j);
                neutral_stable_g01(end+1) = g01_root_vals(j);
                neutral_stable_g11(end+1) = g11_root_vals(j);
                neutral_stable_s(end+1) = s_val_range(i);
            end
        end
    end
end

neutral_stable_q = neutral_stable_g11 + neutral_stable_g01;
selected_stable_q = selected_stable_g11 + selected_stable_g01;
unstable_q = unstable_g11 + unstable_g01;