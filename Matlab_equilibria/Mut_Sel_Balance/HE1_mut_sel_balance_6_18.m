% for the HE1 case under selection-mutation balance

syms g00 g01 g10 g11 q p pa pb qa qb s h1 h2 h3 mu H D G F

assume(q>=0 & q<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(D>=0 & D<=1/4);


% equations to add selection
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

% simplifying equation of linkage disequilibrium
% simplifies to g01 = -q^2 + q - D
linkage_disequilibrium = D == g00*g11 - g01*g10;
ld_2 = subs(linkage_disequilibrium, g11, 1-(g00+g01+g10));
ld_3 = subs(ld_2, g00, q-g10);
ld_4 = subs(ld_3, g10, g01);
ld_5 = isolate(ld_4, g01);

%heterozygosity = H == g01 + g10 + D;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 - g00 == 0;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu) - g01 == 0;
mut_g10 = sel_g00*mu*(1-mu) + sel_g10*(1-mu) - g10 == 0;
mut_g11 = sel_g00*mu^2 + sel_g01*mu + sel_g10*mu + sel_g11 - g11 == 0;


mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g11 from the equation by replacing it with 1-(g00+g01+g10)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g11, 1-(g00+g01+g10));
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g11, (p^2 + D));
    
    % removes g00 from the equation using the relationship that q = g10 + g00 (implicitly that g10 = g01)
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g00, q-g10);
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g00, (q^2 + D));
    
    % removes g10 from the equation using g01 = g10
    mut_eqn_set(i) = subs(mut_eqn_set(i), g10, g01);
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g10, g01);
    
    % removes g01 from the equation using above expression (ld_5) for D
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g01, -q^2 + q - D);
    %mut_eqn_set(i) = subs(mut_eqn_set(i), g01, (p*q + D));
end


%Y = solve(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), s, 'ReturnConditions', true)

%Y_2 = vpasolve([mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4)], q, .01)
%eqn5 = g00+g01+g10+g11 == 1; - used in simplifications above using subs()
%eqn6 = g00+g01 == qa; - used in simplifications above using subs()
%eqn7 = g00+g10 == qb; - used in simplifications above using subs()
%eqn8 = g10 == g01; - used in simplifications above using subs()
