% for the HE2 case

syms g00 g01 g10 g11 pa pb qa qb

eqn1 = g00^2+(9/8)*g00*g01+(9/8)*g00*g10+(1/4)*g01^2+(1/4)*g10^2+(1/2)*g01*g10+(1/2)*g00*g11+(1/8)*g01*g11+(1/8)*g10*g11 == g00;

eqn2 = (3/8)*g00*g01+(3/8)*g00*g10+(1/4)*g01^2+(1/4)*g10^2+(1/2)*g01*g10+(1/2)*g00*g11+(3/8)*g01*g11+(3/8)*g10*g11 == g10;

eqn3 = (3/8)*g00*g01+(3/8)*g00*g10+(1/4)*g01^2+(1/4)*g10^2+(1/2)*g01*g10+(1/2)*g00*g11+(3/8)*g01*g11+(3/8)*g10*g11 == g01;

eqn4 = (1/8)*g00*g01+(1/8)*g00*g10+(1/4)*g01^2+(1/4)*g10^2+(1/2)*g01*g10+(1/2)*g00*g11+(9/8)*g01*g11+(9/8)*g10*g11+g11^2 == g11;

eqn5 = g00+g01+g10+g11 == 1;
eqn6 = g11+g10 == pa;
eqn7 = g11+g01 == pb;
eqn8 = qa == 1-pa;
eqn9 = qb == 1-pb;
eqn10 = pa == pb; % equivalent to g01 == g10 and also implies qa == qb (this is the key simplification)

% below is equivalent to stating that the sum of the genotype frequencies
% must equal one... which is an unnecessary equation, as it turns out
eqn11 = g00^2 + (2*g00*g01 + 2*g00*g10) + (g10^2 + g01^2 + 2*(g00*g11 + g01*g10)) + (2*g01*g11 + 2*g10*g11) + g11^2 == 1; 

%eqn_system = [eqn1 eqn2 eqn3 eqn4 eqn5];

S = solve(eqn1, eqn2, eqn3, eqn5, eqn6, eqn7, 'ReturnConditions', true)


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