% for the HE2 case

syms g00 g01 g10 g11 p q pa pb qa qb beta gamma

assume(gamma>=0 & gamma<=1);
assume(beta>=0 & beta<=1);
assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);
assume(q>=0 & q<=1);


meiosis_g00 = g00^2 + (.5+(beta+beta*gamma)/32)*2*g01*g00 + (.5+(beta+beta*gamma)/32)*2*g10*g00 + ((3*beta+beta*gamma)/16)*g01^2 + ((3*beta+beta*gamma)/16)*g10^2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + ((beta+beta*gamma)/32)*2*g01*g11 + ((beta+beta*gamma)/32)*2*g10*g11 == g00;
meiosis_g10 = ((3*beta+3*beta*gamma)/32)*2*g01*g00 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g00 + ((beta+3*beta*gamma)/16)*g01^2 + (1-(7*beta+5*beta*gamma)/16)*g10^2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + ((3*beta+3*beta*gamma)/32)*2*g01*g11 + (.5-(5*beta+5*beta*gamma)/32)*2*g10*g11 == g10;
meiosis_g01 = (.5-(5*beta+5*beta*gamma)/32)*2*g01*g00 + ((3*beta+3*beta*gamma)/32)*2*g10*g00 + (1-(7*beta+5*beta*gamma)/16)*g01^2 + ((beta+3*beta*gamma)/16)*g10^2 + (.25+(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + (.5-(5*beta+5*beta*gamma)/32)*2*g01*g11 + ((3*beta+3*beta*gamma)/32)*2*g10*g11 == g01;
meiosis_g11 = ((beta+beta*gamma)/32)*2*g01*g00 + ((beta+beta*gamma)/32)*2*g10*g00 + ((3*beta+beta*gamma)/16)*g01^2 + ((3*beta+beta*gamma)/16)*g10^2 + (.25-(beta-beta*gamma)/32)*(2*g00*g11+2*g01*g10) + (.5+(beta+beta*gamma)/32)*2*g01*g11 + (.5+(beta+beta*gamma)/32)*2*g10*g11 + g11^2 == g11;

eqn5 = g00+g01+g10+g11 == 1;
eqn6 = g11+g10 == q;
eqn7 = g11+g01 == q;
eqn10 = g01 == g10;
eqn8 = pa == 1-qa;
eqn9 = pb == 1-qb;
eqn12 = qa == qb; % equivalent to g01 == g10 and also implies qa == qb (this is the key simplification)

% below is equivalent to stating that the sum of the genotype frequencies
% must equal one... which is an unnecessary equation, as it turns out
eqn13 = g00^2 + (2*g00*g01 + 2*g00*g10) + (g10^2 + g01^2 + 2*(g00*g11 + g01*g10)) + (2*g01*g11 + 2*g10*g11) + g11^2 == 1; 

%eqn_system = [eqn1 eqn2 eqn3 eqn4 eqn5];

S = solve(meiosis_g00, meiosis_g10, meiosis_g01, meiosis_g11, eqn5, eqn6, eqn7, eqn10, g00, g01, g10, g11, 'ReturnConditions', true)

G00 = factor(expand(S.g00^2));
G01 = factor(expand(2*S.g00*S.g01));
G10 = factor(expand(2*S.g00*S.g10));
G02 = factor(expand(S.g01^2));
G20 = factor(expand(S.g10^2));
G11 = factor(expand(2*(S.g11*S.g00 + S.g01*S.g10)));
G21 = factor(expand(2*S.g11*S.g10));
G12 = factor(expand(2*S.g11*S.g01));
G22 = factor(expand(S.g11^2));

G0 = G00
G1 = factor(expand(2*S.g00*S.g01 + 2*S.g00*S.g10))
G2 = factor(expand(S.g01^2 + S.g10^2 + 2*(S.g11*S.g00 + S.g01*S.g10)))
G3 = factor(expand(2*S.g11*S.g01 + 2*S.g11*S.g10))
G4 = G22
