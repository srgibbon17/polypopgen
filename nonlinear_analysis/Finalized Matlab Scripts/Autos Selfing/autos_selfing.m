% for autos, classification of fixed points using linear stability
% analysis, the Jacobian matrix, and eigenvectors

iterations = 5; % number of steps or number of data points to generate

%s_val_range = logspace(-8, -5, iterations); % starting s value
s_val_range = 1e-7;

mu_val = 2e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
a_val = 0; % constant value of alpha (double reduction rate)
sigma_val = 0; % selfing rate

h1_val = 0; % h1 dominance coefficient value, constant
h2_val = 0; % h2 dominance coefficient value, constant
h3_val = 0; % h3 dominance coefficient value, constant

syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu sigma

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
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2);
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2);

G0_t_1 = mut_g0^2;
G1_t_1 = 2*mut_g0*mut_g1;
G2_t_1 = 2*mut_g0*mut_g2 + mut_g1^2;
G3_t_1 = 2*mut_g1*mut_g2;
G4_t_1 = mut_g2^2;

delta_G0 = G0_t_1 - G0;
delta_G1 = G1_t_1 - G1;
delta_G2 = G2_t_1 - G2;
delta_G3 = G3_t_1 - G3;
delta_G4 = G4_t_1 - G4;

genotypes = [delta_G0, delta_G1, delta_G2, delta_G3, delta_G4];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(genotypes)
    genotypes(i) = subs(genotypes(i), G2, (1-G0-G1-G3-G4));
end

initial_guess = [0.014588449116675992046621873242267, 0.1095524452929494429669727926448, 0.38612504070035016963499313612225, 0.18122616803745835201520938256452];

for i = 1:length(s_val_range)
    [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = root_solns(genotypes(1), genotypes(2), genotypes(4), genotypes(5), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, G0, G1, G3, G4, initial_guess)
end

% for i = 1:length(s_val_range)
% 
%     %solves for the fixed points of the system
%     [g0_root_vals, g1_root_vals] = root_solns(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);
% 
%     %evaluating the jacobian and stability of each fixed point
%     [fixed_pt_stabilities] = linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals); 
% 
%     for j = 1:length(fixed_pt_stabilities)
% 
%         if fixed_pt_stabilities(j) == 0
%             unstable_g0(end+1) = g0_root_vals(j);
%             unstable_g1(end+1) = g1_root_vals(j);
%             unstable_s(end+1) = s_val_range(i);
% 
%         elseif fixed_pt_stabilities(j) == 1
%             if g0_root_vals(j) > .3333
%                 selected_stable_g0(end+1) = g0_root_vals(j);
%                 selected_stable_g1(end+1) = g1_root_vals(j);
%                 selected_stable_s(end+1) = s_val_range(i);
%             else
%                 neutral_stable_g0(end+1) = g0_root_vals(j);
%                 neutral_stable_g1(end+1) = g1_root_vals(j);
%                 neutral_stable_s(end+1) = s_val_range(i);
%             end
%         end
%     end
% end

% figure
% 
% plot(neutral_stable_s, neutral_stable_g0+.5*neutral_stable_g1)
% hold on
% plot(selected_stable_s, selected_stable_g0+.5*selected_stable_g1)
% plot(unstable_s, unstable_g0+.5*unstable_g1, 'LineStyle','--')
% 
% xscale log


%%% FUNCTIONS %%%

function [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = root_solns(mut_G0_eqn, mut_G1_eqn, mut_G3_eqn, mut_G4_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, G0, G1, G3, G4, initial_guess)

    %function which uses vpasolve to find the fixed points/roots of the system

    G0_eqn = subs(mut_G0_eqn, mu, mu_val);
    G0_eqn = subs(G0_eqn, nu, nu_val);
    G0_eqn = subs(G0_eqn, s, s_val);
    G0_eqn = subs(G0_eqn, h1, h1_val);
    G0_eqn = subs(G0_eqn, h2, h2_val);
    G0_eqn = subs(G0_eqn, h3, h3_val);
    G0_eqn = subs(G0_eqn, a, a_val);

    G1_eqn = subs(mut_G1_eqn, mu, mu_val);
    G1_eqn = subs(G1_eqn, nu, nu_val);
    G1_eqn = subs(G1_eqn, s, s_val);
    G1_eqn = subs(G1_eqn, h1, h1_val);
    G1_eqn = subs(G1_eqn, h2, h2_val);
    G1_eqn = subs(G1_eqn, h3, h3_val);
    G1_eqn = subs(G1_eqn, a, a_val);

    G3_eqn = subs(mut_G3_eqn, mu, mu_val);
    G3_eqn = subs(G3_eqn, nu, nu_val);
    G3_eqn = subs(G3_eqn, s, s_val);
    G3_eqn = subs(G3_eqn, h1, h1_val);
    G3_eqn = subs(G3_eqn, h2, h2_val);
    G3_eqn = subs(G3_eqn, h3, h3_val);
    G3_eqn = subs(G3_eqn, a, a_val);

    G4_eqn = subs(mut_G4_eqn, mu, mu_val);
    G4_eqn = subs(G4_eqn, nu, nu_val);
    G4_eqn = subs(G4_eqn, s, s_val);
    G4_eqn = subs(G4_eqn, h1, h1_val);
    G4_eqn = subs(G4_eqn, h2, h2_val);
    G4_eqn = subs(G4_eqn, h3, h3_val);
    G4_eqn = subs(G4_eqn, a, a_val);


    [G0_root_vals, G1_root_vals, G3_root_vals, G4_root_vals] = vpasolve([G0_eqn, G1_eqn, G3_eqn, G4_eqn], [G0, G1, G3, G4], initial_guess);
end


function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_val, g1, g1_root_val)
    
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

function [fixed_pt_stabilities] = linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals, g1, g1_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, det(J) and tr(J), classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g0_root_vals;

    jacobian_eval = [0, 0; 0, 0];

    for i = 1:length(g0_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = pd_evaluation(jacobian_matrix(j, k), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_root_vals(i), g1, g1_root_vals(i)); 
            end
        end

        %calulating the trace and determinant of the evaluated jacobian
        trace_jac = trace(jacobian_eval);
        det_jac = det(jacobian_eval);

        [eigenvectors, eigenvalues] = eig(jacobian_eval);

        disp(eigenvalues)
        disp(eigenvectors)

        %eigenvals = [abs(eigenvalues(1,1)), abs(eigenvalues(2,2))];

        %disp(max(eigenvals)/min(eigenvals))

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
