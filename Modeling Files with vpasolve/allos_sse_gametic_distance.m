syms a b c d qa qb

g_t = [1-b-c-d, b; c, d];
g_t_1 = [a^2+a*b+a*c+(1/2)*a*d+(1/2)*b*c, b^2+a*b+b*d+(1/2)*a*d+(1/2)*b*c;  c^2+a*c+c*d+(1/2)*a*d+(1/2)*b*c, d^2+b*d+c*d+(1/2)*a*d+(1/2)*b*c];

g_star = [(1-qa)*(1-qb), (1-qa)*qb; (1-qb)*qa, qa*qb];

for i = 1:length(g_t_1)
    for j = 1:length(g_t_1)
    g_t_1(i,j) = subs(g_t_1(i,j), a, 1-b-c-d);
    g_t_1(i,j) = simplify(expand(g_t_1(i,j)));
    end
end

for i = 1:length(g_star)
    for j = 1:length(g_star)
    g_star(i,j) = subs(g_star(i,j), qa, c+d);
    g_star(i,j) = subs(g_star(i,j), qb, b+d);
    g_star(i,j) = simplify(expand(g_star(i,j)));
    end
end

g_diff_t_1_star = g_t_1 - g_star;

g_diff_t_1_star_squared = g_diff_t_1_star.*g_diff_t_1_star;

g_t_1_star_sum = 0;

for i = 1:length(g_diff_t_1_star_squared)
    for j = 1:length(g_diff_t_1_star_squared)
        g_t_1_star_sum = g_t_1_star_sum + g_diff_t_1_star_squared(i,j);
    end
end

simplify(expand(g_t_1_star_sum));

g_diff_t_star = g_t - g_star;

g_diff_t_star_squared = g_diff_t_star.*g_diff_t_star;

g_t_star_sum = 0;

for i = 1:length(g_diff_t_star_squared)
    for j = 1:length(g_diff_t_star_squared)
        g_t_star_sum = g_t_star_sum + g_diff_t_star_squared(i,j);
    end
end

simplify(expand(g_t_star_sum))

variable = simplify(expand(g_t_1_star_sum/g_t_star_sum));

