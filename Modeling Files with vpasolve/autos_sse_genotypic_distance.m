syms a b c d e q

G_t = [1-b-c-d-e, b, c, d, e];

g_t_1 = [a+(1/2)*b+(1/6)*c, (1/2)*b+(2/3)*c+(1/2)*d, (1/6)*c+(1/2)*d+e];

G_t_1 = [(a+(1/2)*b+(1/6)*c)^2, 2*(a+(1/2)*b+(1/6)*c)*((1/2)*b+(2/3)*c+(1/2)*d), ((1/2)*b+(2/3)*c+(1/2)*d)^2 + 2*(a+(1/2)*b+(1/6)*c)*((1/6)*c+(1/2)*d+e), 2*((1/2)*b+(2/3)*c+(1/2)*d)*((1/6)*c+(1/2)*d+e), ((1/6)*c+(1/2)*d+e)^2 ];

G_star = [(1-q)^4, 4*q*(1-q)^3, 6*q^2*(1-q)^2, 4*q^3*(1-q), q^4];

for i = 1:length(G_t_1)
    G_t_1(i) = subs(G_t_1(i), a, 1-b-c-d-e);
    %G_t_1(i) = simplify(expand(G_t_1(i)));
end

for i = 1:length(G_star)
    G_star(i) = subs(G_star(i), q, (1/4)*b+(1/2)*c+(3/4)*d+e);
    %G_star(i) = simplify(expand(G_star(i)));
end

G_diff_t_1_star = G_t_1 - G_star;

G_diff_t_1_star_squared = G_diff_t_1_star.*G_diff_t_1_star;

G_t_1_star_sum = 0;

for i = 1:length(G_diff_t_1_star_squared)
    G_t_1_star_sum = G_t_1_star_sum + G_diff_t_1_star_squared(i);
end

G_diff_t_star = G_t - G_star;

G_diff_t_star_squared = G_diff_t_star.*G_diff_t_star;

G_t_star_sum = 0;

for i = 1:length(G_diff_t_star_squared)
    G_t_star_sum = G_t_star_sum + G_diff_t_star_squared(i);
end

simplify(expand(G_t_1_star_sum))

simplify(expand(G_t_star_sum))
