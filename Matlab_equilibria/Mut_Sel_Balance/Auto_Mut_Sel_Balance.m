% for the autopolyploids

syms alpha s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 w0 w1 w2 w3 w4 wbar mu D

wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);

w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

sel_meiosis_g0 = G0*w0+(1/2 + alpha/4)*G1*w1 + (1/6 + alpha/3)*G2*w2 + (alpha/4)*G3*w3;
sel_meiosis_g1 = (1/2 - alpha/2)*G1*w1 + (2/3 - 2*alpha/3)*G2*w2 + (1/2 - alpha/2)*G3*w3;
sel_meiosis_g2 = (alpha/4)*G1*w1 + (1/6 + alpha/3)*G2*w2 + (1/2 + alpha/4)*G3*w3 + G4*w4;

mut_g0 = sel_meiosis_g0*(1-mu)^2 - g0 == 0;
mut_g1 = 2*sel_meiosis_g0*(1-mu)*mu + sel_meiosis_g1*(1-mu) - g1 == 0;
mut_g2 = sel_meiosis_g0*mu^2 + sel_meiosis_g1*mu + sel_meiosis_g2 - g2 == 0;

eqn_set = [mut_g0, mut_g1, mut_g2];

for i = 1:length(eqn_set)
    eqn_set(i) = subs(eqn_set(i), G0, g0^2);
    eqn_set(i) = subs(eqn_set(i), G1, 2*g0*g1);
    eqn_set(i) = subs(eqn_set(i), G2, (2*g0*g2 + g1^2));
    eqn_set(i) = subs(eqn_set(i), G3, 2*g1*g2);
    eqn_set(i) = subs(eqn_set(i), G4, g2^2);
    eqn_set(i) = subs(eqn_set(i), g2, (1-2*g1-g0));
    eqn_set(i) = subs(eqn_set(i), g1, (q-g0));
    % using a measure for D
    %eqn_set(i) = subs(eqn_set(i), g0, (q + ((16*q^2)/3 - (16*q)/3 + 4*D + 16/9)^(1/2)/2 - 2/3));
    % using a measure for H
    eqn_set(i) = subs(eqn_set(i), g0, (q + ((16*q^2)/3 - (16*q)/3 + 4*H + 1/9)^(1/2)/2 - 1/6));
end

Y = solve(eqn_set(1), eqn_set(2), eqn_set(3), mu, 'ReturnConditions', true)

% g0 + 2g1 + g2 = 0
% q = g0 + g1
% D = g0g2 - g1^2