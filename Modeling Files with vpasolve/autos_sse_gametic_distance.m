syms a b c q alpha

g_t = [1-b-c, b, c];
g_t_1 = [a^2+a*b+(1/6)*b^2+(1/3)*a*c, a*b+(2/3)*b^2+(4/3)*a*c+b*c,  (1/6)*b^2+(1/3)*a*c+b*c+c^2];

g_star = [(1-q)^2, 2*(1-q)*q, q^2];

for i = 1:length(g_t_1)
    g_t_1(i) = subs(g_t_1(i), a, 1-b-c);
    g_t_1(i) = simplify(expand(g_t_1(i)));
end

for i = 1:length(g_star)
    g_star(i) = subs(g_star(i), q, (1/2)*b+c);
    g_star(i) = simplify(expand(g_star(i)));
end

g_diff_t_1_star = g_t_1 - g_star;

g_diff_t_1_star_squared = g_diff_t_1_star.*g_diff_t_1_star;

g_t_1_star_sum = 0;

for i = 1:length(g_diff_t_1_star_squared)
    g_t_1_star_sum = g_t_1_star_sum + g_diff_t_1_star_squared(i);
end

simplify(expand(g_t_1_star_sum));

g_diff_t_star = g_t - g_star;

g_diff_t_star_squared = g_diff_t_star.*g_diff_t_star;

g_t_star_sum = 0;

for i = 1:length(g_diff_t_star_squared)
    g_t_star_sum = g_t_star_sum + g_diff_t_star_squared(i);
end

simplify(expand(g_t_star_sum));

simplify(expand(g_t_1_star_sum/g_t_star_sum))

