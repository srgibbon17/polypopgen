% plotting script to generate a figure with 3x3 subplots which are each a
% bifurcation diagram with 3 lines

iterations = 5; 

grid_row = 3;
grid_col = 5;

s_val_range = logspace(-9, -4, iterations); % set of selection coefficients

mu_val = 2e-8; % forward mutation rate
nu_val = 1e-9; % backward mutation rate
a_val = 0; % double reduction rate

a_val_1 = 0;
a_val_2 = 1/12;
a_val_3 = 1/6;

%defining variables so HPC Matlab won't crash (fingers crossed)

jacobian_entry = 0;
jacobian_matrix = 0;
h_val = 0;
h1_val = 0;
h2_val = 0;
h3_val = 0;
s_val = 0;
mut_g00_eqn = 0;
mut_g01_eqn = 0;
mut_g11_eqn = 0;
grid_pos = 1;
mut_eqn_g0 = 0;
mut_g0_eqn = 0;
mut_g1_eqn = 0;
mut_g00_eqn = 0;
mut_g01_eqn = 0;
mut_g11_eqn = 0;
g0_roots = 0;
g0_derivative = 0;
g0_root_val = 0;
g0_root_vals = 0;
g1_root_val = 0;
g1_root_vals = 0;
g00_root_val = 0;
g00_root_vals = 0;
g01_root_val = 0;
g01_root_vals = 0;
g11_root_val = 0;
g11_root_vals = 0;

syms mu nu h1 h2 h3 h g0 g1 g2 g00 g01 g10 g11 a

% dominance coefficients from Kacser and Burns
function [h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val)
    
    % evaluates the three dominance coefficients for various values of k
    % assuming d = 1

    h1_val = (1/3) / (1/3 + (k_val/(1-k_val)));

    h2_val = 1 / (1 + (k_val/(1-k_val)));

    h3_val = 3 / (3 + (k_val/(1-k_val)));

end


% Subfunctions for data generation
%allos:
function [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, mu, mu_val);
    g00_eqn = subs(g00_eqn, nu, nu_val);
    g00_eqn = subs(g00_eqn, s, s_val);
    g00_eqn = subs(g00_eqn, h1, h1_val);
    g00_eqn = subs(g00_eqn, h2, h2_val);
    g00_eqn = subs(g00_eqn, h3, h3_val);

    g01_eqn = subs(mut_g01_eqn, mu, mu_val);
    g01_eqn = subs(g01_eqn, nu, nu_val);
    g01_eqn = subs(g01_eqn, s, s_val);
    g01_eqn = subs(g01_eqn, h1, h1_val);
    g01_eqn = subs(g01_eqn, h2, h2_val);
    g01_eqn = subs(g01_eqn, h3, h3_val);

    g11_eqn = subs(mut_g11_eqn, mu, mu_val);
    g11_eqn = subs(g11_eqn, nu, nu_val);
    g11_eqn = subs(g11_eqn, s, s_val);
    g11_eqn = subs(g11_eqn, h1, h1_val);
    g11_eqn = subs(g11_eqn, h2, h2_val);
    g11_eqn = subs(g11_eqn, h3, h3_val);


    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end

function [pd_value] = allo_pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_val, g01, g01_root_val, g11, g11_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, s, s_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, g00, g00_root_val);
    pd_value = subs(pd_value, g01, g01_root_val);
    pd_value = subs(pd_value, g11, g11_root_val);
end

function [fixed_pt_stabilities] = allo_linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g00_root_vals;

    jacobian_eval = [0, 0, 0; 
                     0, 0, 0;
                     0, 0, 0];

    for i = 1:length(g00_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = allo_pd_evaluation(jacobian_matrix(j, k), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals(i), g01, g01_root_vals(i), g11, g11_root_vals(i)); 
            end
        end

        %calulating the eigenvalues of the evaluated jacobian
        [eigenvectors, eigenvalues] = eig(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if eigenvalues(1,1) < 0 && eigenvalues(2,2) < 0 && eigenvalues(3,3) < 0
            fixed_pt_stabilities(i) = 1; %0 for stable node
        elseif eigenvalues(1,1) > 0 && eigenvalues(2,2) > 0 && eigenvalues(3,3) > 0
            disp('Error. Unexpected stability type from linear stability analysis.')
        else
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        end
    end
end

%autos:
function [g0_root_vals, g1_root_vals] = auto_root_solns(mut_g0_eqn, mut_g1_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1)

    %function which uses vpasolve to find the fixed points/roots of the system

    g0_eqn = subs(mut_g0_eqn, mu, mu_val);
    g0_eqn = subs(g0_eqn, nu, nu_val);
    g0_eqn = subs(g0_eqn, s, s_val);
    g0_eqn = subs(g0_eqn, h1, h1_val);
    g0_eqn = subs(g0_eqn, h2, h2_val);
    g0_eqn = subs(g0_eqn, h3, h3_val);
    g0_eqn = subs(g0_eqn, a, a_val);

    g1_eqn = subs(mut_g1_eqn, mu, mu_val);
    g1_eqn = subs(g1_eqn, nu, nu_val);
    g1_eqn = subs(g1_eqn, s, s_val);
    g1_eqn = subs(g1_eqn, h1, h1_val);
    g1_eqn = subs(g1_eqn, h2, h2_val);
    g1_eqn = subs(g1_eqn, h3, h3_val);
    g1_eqn = subs(g1_eqn, a, a_val);


    [g0_root_vals, g1_root_vals] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);
end

function [pd_value] = auto_pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_val, g1, g1_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, s, s_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, a, a_val);
    pd_value = subs(pd_value, g0, g0_root_val);
    pd_value = subs(pd_value, g1, g1_root_val);
end

function [fixed_pt_stabilities] = auto_linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, det(J) and tr(J), classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g0_root_vals;

    jacobian_eval = [0, 0; 0, 0];

    for i = 1:length(g0_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = auto_pd_evaluation(jacobian_matrix(j, k), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals(i), g1, g1_root_vals(i)); 
            end
        end

        %calulating the trace and determinant of the evaluated jacobian
        trace_jac = trace(jacobian_eval);
        det_jac = det(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if det_jac < 0
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        elseif det_jac > 0
            if trace_jac < 0 && trace_jac^2 - 4*det_jac > 0
                fixed_pt_stabilities(i) = 1; %1 for stable node
            else
                disp('Error. Unexpected stability type from linear stability analysis.')
            end
        else
            disp('Error. Unexpected stability type from linear stability analysis.')
        end
    end
end
%diploids:
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

% Functions for data generation and sorting: 

% auto data ---------------------------------------------------------------
function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val)
%Generates autopolyploid bifurcation plotting data
%   for autos, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       6 data arrays in the following order:
%           neutral_stable (for q and s, respectively)
%           selected_stable (for q and s, respectively)
%           unstable (for q and s, respectively)

syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=-1 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(a>=0 & a<=1/6);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);
assume(G3>=0 & G3<=1);
assume(G4>=0 & G4<=1);

% equations to parameterize relative fitnesses
wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g0 = G0*w0+(1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3;
sel_g1 = (1/2 - a/2)*G1*w1 + (2/3 - 2*a/3)*G2*w2 + (1/2 - a/2)*G3*w3;
sel_g2 = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4;

% equations for mutation
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0;
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu - g1;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2;

mut_exp_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), G0, g0^2);
    mut_exp_set(i) = subs(mut_exp_set(i), G1, 2*g0*g1);
    mut_exp_set(i) = subs(mut_exp_set(i), G2, (2*g0*g2 + g1^2));
    mut_exp_set(i) = subs(mut_exp_set(i), G3, 2*g1*g2);
    mut_exp_set(i) = subs(mut_exp_set(i), G4, g2^2);
    mut_exp_set(i) = subs(mut_exp_set(i), g2, (1-g1-g0));
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_exp_set(1), g0), diff(mut_exp_set(1), g1); diff(mut_exp_set(2), g0), diff(mut_exp_set(2), g1)];

neutral_stable_g0 = [];
neutral_stable_g1 = [];
neutral_stable_s = [];

selected_stable_g0 = [];
selected_stable_g1 = [];
selected_stable_s = [];

unstable_g0 = [];
unstable_g1 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g0_root_vals, g1_root_vals] = auto_root_solns(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = auto_linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals); 

    for j = 1:length(fixed_pt_stabilities)

        if fixed_pt_stabilities(j) == 0
            unstable_g0(end+1) = g0_root_vals(j);
            unstable_g1(end+1) = g1_root_vals(j);
            unstable_s(end+1) = s_val_range(i);

        elseif fixed_pt_stabilities(j) == 1
            if g0_root_vals(j) > .3333
                selected_stable_g0(end+1) = g0_root_vals(j);
                selected_stable_g1(end+1) = g1_root_vals(j);
                selected_stable_s(end+1) = s_val_range(i);
            else
                neutral_stable_g0(end+1) = g0_root_vals(j);
                neutral_stable_g1(end+1) = g1_root_vals(j);
                neutral_stable_s(end+1) = s_val_range(i);
            end
        end
    end
end

neutral_stable_q = ones(1, length(neutral_stable_g0)) - neutral_stable_g0 - .5*neutral_stable_g1;
selected_stable_q = ones(1, length(selected_stable_g0)) - selected_stable_g0 - .5*selected_stable_g1;
unstable_q = ones(1, length(unstable_g0)) - unstable_g0 - .5*unstable_g1;


end

% diploid data ------------------------------------------------------------
function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val)
%Generates diploid bifurcation plotting data
%   for diploids, a nonlinear model and analysis
%   returns:
%       diploid_data:
%           a 6 element array of data arrays in the following order:
%               neutral_stable (for q and s, respectively)
%               selected_stable (for q and s, respectively)
%               unstable (for q and s, respectively)

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
            if g0_roots(j) > .3333 && real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                selected_stable_g0(end+1) = g0_roots(j);
                selected_stable_s(end+1) = s_val_range(i);
            elseif real(g0_roots(j)) <= 1 && real(g0_roots(j)) >= 0
                neutral_stable_g0(end+1) = g0_roots(j);
                neutral_stable_s(end+1) = s_val_range(i);
            end
        end
    end
end

neutral_stable_q = ones(1, length(neutral_stable_g0)) - neutral_stable_g0;
selected_stable_q = ones(1, length(selected_stable_g0)) - selected_stable_g0;
unstable_q = ones(1, length(unstable_g0)) - unstable_g0;

end

%allo data (incl. general subfunctions) -----------------------------------
function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val)
%Generates 0HE bifurcation plotting data
%   for allos with 0 HEs, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       HE0_data:
%           a 6 element array of data arrays in the following order:
%               neutral_stable (for q and s, respectively)
%               selected_stable (for q and s, respectively)
%               unstable (for q and s, respectively)

syms g00 g01 g10 g11 s h1 h2 h3 mu nu

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
sel_g00 = g00^2*w0+g00*g01*w1+g00*g10*w1+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2;
sel_g10 = g00*g10*w1+g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g10*g11*w3;
sel_g01 = g00*g01*w1+g01^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3;
sel_g11 = (1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3+g10*g11*w3+g11^2*w4;

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
    [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = allo_linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);

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

end

function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = HE1_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val)
%Generates 0HE bifurcation plotting data
%   for allos with 0 HEs, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       HE0_data:
%           a 6 element array of data arrays in the following order:
%               neutral_stable (for q and s, respectively)
%               selected_stable (for q and s, respectively)
%               unstable (for q and s, respectively)

syms g00 g01 g10 g11 s h1 h2 h3 mu nu

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
sel_g00 = g00^2*w0+(17/16)*g00*g01*w1+(17/16)*g00*g10*w1+(3/16)*g01^2*w2+(3/16)*g10^2*w2+(7/16)*g01*g10*w2+(7/16)*g00*g11*w2+(1/16)*g01*g11*w3+(1/16)*g10*g11*w3;
sel_g10 = (3/16)*g00*g01*w1+(11/16)*g00*g10*w1+(1/16)*g01^2*w2+(9/16)*g10^2*w2+(9/16)*g01*g10*w2+(9/16)*g00*g11*w2+(3/16)*g01*g11*w3+(11/16)*g10*g11*w3;
sel_g01 = (11/16)*g00*g01*w1+(3/16)*g00*g10*w1+(9/16)*g01^2*w2+(1/16)*g10^2*w2+(9/16)*g01*g10*w2+(9/16)*g00*g11*w2+(11/16)*g01*g11*w3+(3/16)*g10*g11*w3;
sel_g11 = (1/16)*g00*g01*w1+(1/16)*g00*g10*w1+(3/16)*g01^2*w2+(3/16)*g10^2*w2+(7/16)*g01*g10*w2+(7/16)*g00*g11*w2+(17/16)*g01*g11*w3+(17/16)*g10*g11*w3+g11^2*w4;

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
    [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = allo_linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);

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

end

function [neutral_stable_q, neutral_stable_s, selected_stable_q, selected_stable_s, unstable_q, unstable_s] = HE2_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val)
%Generates 0HE bifurcation plotting data
%   for allos with 0 HEs, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors
%   returns:
%       HE0_data:
%           a 6 element array of data arrays in the following order:
%               neutral_stable (for q and s, respectively)
%               selected_stable (for q and s, respectively)
%               unstable (for q and s, respectively)

syms g00 g01 g10 g11 s h1 h2 h3 mu nu

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
sel_g00 = g00^2*w0+(9/8)*g00*g01*w1+(9/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(1/8)*g01*g11*w3+(1/8)*g10*g11*w3;
sel_g10 = (3/8)*g00*g01*w1+(3/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(3/8)*g01*g11*w3+(3/8)*g10*g11*w3;
sel_g01 = (3/8)*g00*g01*w1+(3/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(3/8)*g01*g11*w3+(3/8)*g10*g11*w3;
sel_g11 = (1/8)*g00*g01*w1+(1/8)*g00*g10*w1+(1/4)*g01^2*w2+(1/4)*g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+(9/8)*g01*g11*w3+(9/8)*g10*g11*w3+g11^2*w4;

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
    [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = allo_linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals);

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

end

%plotting functions -------------------------------------------------------

function ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [autos_q1, autos_s1, autos_q2, autos_s2, autos_q3, autos_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val);

    [diploids_q1, diploids_s1, diploids_q2, diploids_s2, diploids_q3, diploids_s3] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);

    if h_val <= 2/3

        HE0_q_values = cat(2, HE0_q1, HE0_q2);
        HE0_s_values = cat(2, HE0_s1, HE0_s2);

        autos_q_values = cat(2, autos_q1, autos_q2);
        autos_s_values = cat(2, autos_s1, autos_s2);

        diploids_q_values = cat(2, diploids_q1, diploids_q2);
        diploids_s_values = cat(2, diploids_s1, diploids_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s_values, HE0_q_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(autos_s_values, autos_q_values, 'Color', [0.8500 0.3250 0.0980])
        plot(diploids_s_values, diploids_q_values, 'Color', 'k')

        if h_val == 0
            legend('0 HEs', 'α=0', 'Diploids')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(autos_s1, autos_q1, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s2, autos_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos_s3, autos_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(diploids_s1, diploids_q1, 'Color','k')
        plot(diploids_s2, diploids_q2, 'Color','k')
        plot(diploids_s3, diploids_q3, 'Color','k', 'LineStyle','--')

    end

    xscale log

    if grid_pos == 6 
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 13
        xlabel('s (selection coefficient)')
    end
    
    ylim([0, 1])

    disp(grid_pos)

end

function allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [HE0_q1, HE0_s1, HE0_q2, HE0_s2, HE0_q3, HE0_s3] = HE0_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [HE1_q1, HE1_s1, HE1_q2, HE1_s2, HE1_q3, HE1_s3] = HE1_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    [HE2_q1, HE2_s1, HE2_q2, HE2_s2, HE2_q3, HE2_s3] = HE2_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

    if h_val <= 2/3

        HE0_q_values = cat(2, HE0_q1, HE0_q2);
        HE0_s_values = cat(2, HE0_s1, HE0_s2);

        HE1_q_values = cat(2, HE1_q1, HE1_q2);
        HE1_s_values = cat(2, HE1_s1, HE1_s2);

        HE2_q_values = cat(2, HE2_q1, HE2_q2);
        HE2_s_values = cat(2, HE2_s1, HE2_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s_values, HE0_q_values, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE1_s_values, HE1_q_values, 'Color', [0 0.4470 0.7410 .7])
        plot(HE2_s_values, HE2_q_values, 'Color', [0 0.4470 0.7410 .4])

        if h_val == 0
            legend('0 HEs', '1 HE', '2 HEs')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(HE0_s1, HE0_q1, 'Color', [0 0.4470 0.7410])
        hold on
        plot(HE0_s2, HE0_q2, 'Color', [0 0.4470 0.7410])
        plot(HE0_s3, HE0_q3, 'Color', [0 0.4470 0.7410], 'LineStyle','--')

        plot(HE1_s1, HE1_q1, 'Color', [0 0.4470 0.7410 .7])
        plot(HE1_s2, HE1_q2, 'Color', [0 0.4470 0.7410 .7])
        plot(HE1_s3, HE1_q3, 'Color', [0 0.4470 0.7410 .7], 'LineStyle','--')

        plot(HE2_s1, HE2_q1, 'Color', [0 0.4470 0.7410 .4])
        plot(HE2_s2, HE2_q2, 'Color', [0 0.4470 0.7410 .4])
        plot(HE2_s3, HE2_q3, 'Color', [0 0.4470 0.7410 .4], 'LineStyle','--')

    end
    
    xscale log

    if grid_pos == 6 
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 13 
        xlabel('s (selection coefficient)')
    end

    ylim([0, 1])

    disp(grid_pos)

end

function autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

    [autos1_q1, autos1_s1, autos1_q2, autos1_s2, autos1_q3, autos1_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_1);

    [autos2_q1, autos2_s1, autos2_q2, autos2_s2, autos2_q3, autos2_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_2);
  
    [autos3_q1, autos3_s1, autos3_q2, autos3_s2, autos3_q3, autos3_s3] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val, a_val_3);

    if h_val <= 2/3

        autos_1_q_values = cat(2, autos1_q1, autos1_q2);
        autos_1_s_values = cat(2, autos1_s1, autos1_s2);

        autos_2_q_values = cat(2, autos2_q1, autos2_q2);
        autos_2_s_values = cat(2, autos2_s1, autos2_s2);

        autos_3_q_values = cat(2, autos3_q1, autos3_q2);
        autos_3_s_values = cat(2, autos3_s1, autos3_s2);

        subplot(grid_row, grid_col, grid_pos)

        plot(autos_1_s_values, autos_1_q_values, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(autos_2_s_values, autos_2_q_values, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos_3_s_values, autos_3_q_values, 'Color', [0.8500 0.3250 0.0980 .4])

        if h_val == 0
            legend('α=0', 'α=1/12', 'α=1/6')
        end

    else

        subplot(grid_row, grid_col, grid_pos)

        plot(autos1_s1, autos1_q1, 'Color', [0.8500 0.3250 0.0980])
        hold on
        plot(autos1_s2, autos1_q2, 'Color', [0.8500 0.3250 0.0980])
        plot(autos1_s3, autos1_q3, 'Color', [0.8500 0.3250 0.0980], 'LineStyle','--')

        plot(autos2_s1, autos2_q1, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos2_s2, autos2_q2, 'Color', [0.8500 0.3250 0.0980 .7])
        plot(autos2_s3, autos2_q3, 'Color', [0.8500 0.3250 0.0980 .7], 'LineStyle','--')

        plot(autos3_s1, autos3_q1, 'Color', [0.8500 0.3250 0.0980 .4])
        plot(autos3_s2, autos3_q2, 'Color', [0.8500 0.3250 0.0980 .4])
        plot(autos3_s3, autos3_q3, 'Color', [0.8500 0.3250 0.0980 .4], 'LineStyle','--')

    end
    
    xscale log

    if grid_pos == 6
        ylabel('q (deleterious allele)')
    end

    if grid_pos == 13
        xlabel('s (selection coefficient)')
    end

    ylim([0, 1])

    disp(grid_pos)

end



figure

% fully recessive case ----------------------------------------------------

k_val = 1; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 1;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 6;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 11;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% partially recessive case ------------------------------------------------

k_val = .75; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 2;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 7;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 12;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% additive case -----------------------------------------------------------

k_val = .5; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 3;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 8;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 13;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% partially dominant case -----------------------------------------------------

k_val = .25; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 4;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 9;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 14;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

% fully dominant case -----------------------------------------------------

k_val = 0; % Kacser and Burns parameter

[h1_val, h2_val, h3_val] = sigmoid_dominance_relations(k_val);

h_val = h2_val;

grid_pos = 5;

ploidy_comparison_plot(s_val_range, mu_val, nu_val, a_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 10;

allos_comparison_plot(s_val_range, mu_val, nu_val, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

grid_pos = 15;

autos_comparison_plot(s_val_range, mu_val, nu_val, a_val_1, a_val_2, a_val_3, h_val, h1_val, h2_val, h3_val, grid_row, grid_col, grid_pos)

savefig("15_bifn_diagram_v1.fig")


