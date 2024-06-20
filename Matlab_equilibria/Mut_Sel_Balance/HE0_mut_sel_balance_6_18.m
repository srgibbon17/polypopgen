% for the HE0 case (as a test)

syms g00 g01 g10 g11 pa pb qa qb s h1 h2 h3 mu H D G F

% equations to add selection
wbar = (1-2*s*(h1*(g00*g10+g00*g01)+h2*(g00*g11+g01*g10)+h3*(g01*g11+g10*g11))-s*(h2*(g01^2+g10^2)+g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% selection equations
sel_g00 = g00^2*w0+g00*g01*w1+g00*g10*w1+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2;
sel_g10 = g00*g10*w1+g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g10*g11*w3;
sel_g01 = g00*g01*w1+g01^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3;
sel_g11 = (1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3+g10*g11*w3+g11^2*w4;

sel_eqn_set = [sel_g00, sel_g10, sel_g01, sel_g11];

for i = 1:length(sel_eqn_set)
    % removes g11 from the equation by replacing it with 1-(g00+g01+g10)
    sel_eqn_set(i) = subs(sel_eqn_set(i), g11, 1-(g00+g01+g10));

    % removes g00 from the equation using the relationship that q = g10 + g00 (implicitly that g10 = g01)
    sel_eqn_set(i) = subs(sel_eqn_set(i), g00, q-g10);

    % removes g10 from the equation using g01 = g10
    sel_eqn_set(i) = subs(sel_eqn_set(i), g10, g01);
    
    % removes g01 from the equation using above expression (ld_5) for D
    sel_eqn_set(i) = subs(sel_eqn_set(i), g01, -q^2 + q - D);
end

eqn5 = g00+g01+g10+g11 == 1;


% mutation equations
mut_g00 = sel_g00*(1-mu)^2;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu);
mut_g10 = sel_g00*mu*(1-mu) + sel_g10*(1-mu);
mut_g11 = sel_g00*mu^2 + sel_g01*mu + sel_g10*mu + sel_g11;

eqn6 = g00+g01 == qa;
eqn7 = g00+g10 == qb;
%eqn8 = qa == qb;
%eqn8 = qa == 1-pa;
%eqn9 = qb == 1-pb;
%eqn10 = pa == pb; % equivalent to g01 == g10 and also implies qa == qb

% optional equations to add the assumption of dominance additivity
%eqn16 = h1 == .25;
%eqn17 = h3 == 3*h1;
%eqn18 = h2 == 2*h1;

eqn1 = mut_g00 - g00 == 0;
eqn2 = mut_g01 - g01 == 0;
eqn3 = mut_g10 - g10 == 0;
eqn4 = mut_g11 - g11 == 0;

linkage_disequilibrium = D == g00*g11 - g01*g10;
heterozygosity = H == g01 + g10 + D;
%G = g00 + g11 - D;
%F = g01 - g10; %F-statistic? for difference in allele frequency (q) across subgenomes

%eqn1 = delta_sel_g00 == 0;
%eqn2 = delta_sel_g01 == 0;
%eqn3 = delta_sel_g10 == 0;
%eqn4 = delta_sel_g11 == 0;

Y = solve(eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, 'ReturnConditions', true)
