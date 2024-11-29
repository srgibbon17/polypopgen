function [beta_val, gamma_val] = beta_gamma_calculation(slope, alpha_val)
% Solves exactly for beta and gamma given an alpha value and a slope
% here, alpha corresponds to the double reduction rate and the slope is
% indicates whether synergy or interference is occuring and to what extent
beta = sym('beta');
gamma = sym('gamma');

assume(beta>=0 & beta<=1);
assume(gamma>=0 & gamma<=1);

eqn1 = 1/8*(beta+beta*gamma) == alpha_val;
eqn2 = gamma == slope*beta;

[beta_val, gamma_val] = solve([eqn1, eqn2], [beta, gamma]);

beta_val = double(beta_val);
gamma_val = double(gamma_val);

end