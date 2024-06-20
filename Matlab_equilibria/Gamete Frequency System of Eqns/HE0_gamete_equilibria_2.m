% for the HE0 case (as a test)

syms w x y z m n r s
% w is g00
% x is g01
% y is g10
% z is g11
% m is pa or w+x
% n is pb or w+y

eqn1 = w^2+w*x+w*y+(1/2)*x*y+(1/2)*w*z == w;
eqn2 = w*y+y^2+(1/2)*x*y+(1/2)*w*z+y*z == y;
eqn3 = w*x+x^2+(1/2)*x*y+(1/2)*w*z+x*z == x;
eqn4 = (1/2)*x*y+(1/2)*w*z+x*z+y*z+z^2 == z;
eqn5 = w + x + y + z == 1;

eqn_system = [eqn1 eqn2 eqn3 eqn4 eqn5];

for i = 1:length(eqn_system)
    eqn_system(i) = subs(eqn_system(i), w, m-x);
    eqn_system(i) = subs(eqn_system(i), w, n-y);
    eqn_system(i) = subs(eqn_system(i), z, r-y);
    eqn_system(i) = subs(eqn_system(i), z, s-x);
    eqn_system(i) = subs(eqn_system(i), r, 1-m);
    eqn_system(i) = subs(eqn_system(i), s, 1-n);
    %eqn_system(i) = subs(eqn_system(i), x, y);
end

S = solve(eqn_system(1), eqn_system(2), eqn_system(3), eqn_system(4), eqn_system(5), 'ReturnConditions', true)
