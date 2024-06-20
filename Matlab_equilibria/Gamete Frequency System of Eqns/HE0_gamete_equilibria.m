% for the HE0 case (as a test)

syms g00 g01 g10 g11 pa pb qa qb
% w is g00
% x is g01
% y is g10
% z is g11
% m is pa or w+x
% n is pb or w+y

eqn1 = g00^2+g00*g01+g00*g10+(1/2)*g01*g10+(1/2)*g00*g11 == g00;
eqn2 = g00*g10+g10^2+(1/2)*g01*g10+(1/2)*g00*g11+g10*g11 == g10;
eqn3 = g00*g01+g01^2+(1/2)*g01*g10+(1/2)*g00*g11+g01*g11 == g01;
eqn4 = (1/2)*g01*g10+(1/2)*g00*g11+g01*g11+g10*g11+g11^2 == g11;
eqn5 = g00+g01+g10+g11 == 1;
eqn6 = g11+g10 == pa;
eqn7 = g11+g01 == pb;
eqn8 = qa == 1-pa;
eqn9 = qb == 1-pb;
eqn10 = pa == pb; % equivalent to g01 == g10 and also implies qa == qb
% below is equivalent to stating that the sum of the genotype frequencies must equal one
eqn11 = g00^2 + (2*g00*g01 + 2*g00*g10) + (g10^2 + g01^2 + 2*(g00*g11 + g01*g10)) + (2*g01*g11 + 2*g10*g11) + g11^2 == 1; 

eqn_system = [eqn1 eqn2 eqn3 eqn4 eqn5];

%for i = 1:length(eqn_system)
    %eqn_system(i) = subs(eqn_system(i), x, m-w);
    %eqn_system(i) = subs(eqn_system(i), y, n-w);
    %eqn_system(i) = subs(eqn_system(i), y, r-z);
    %eqn_system(i) = subs(eqn_system(i), x, s-z);
    %eqn_system(i) = subs(eqn_system(i), r, 1-m);
    %eqn_system(i) = subs(eqn_system(i), s, 1-n);
%end

% disp(eqn_system(1))

%S = solve(eqn_system(1), eqn_system(2), eqn_system(3), eqn_system(4), eqn_system(5), 'ReturnConditions', true)

S = solve(eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, 'ReturnConditions', true);

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
