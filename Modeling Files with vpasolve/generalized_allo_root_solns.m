function [g00_root_vals, g01_root_vals, g11_root_vals] = generalized_allo_root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, beta, beta_val, gamma, gamma_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, mu, mu_val);
    g00_eqn = subs(g00_eqn, nu, nu_val);
    g00_eqn = subs(g00_eqn, s, s_val);
    g00_eqn = subs(g00_eqn, h1, h1_val);
    g00_eqn = subs(g00_eqn, h2, h2_val);
    g00_eqn = subs(g00_eqn, h3, h3_val);
    g00_eqn = subs(g00_eqn, beta, beta_val);
    g00_eqn = subs(g00_eqn, gamma, gamma_val);

    g01_eqn = subs(mut_g01_eqn, mu, mu_val);
    g01_eqn = subs(g01_eqn, nu, nu_val);
    g01_eqn = subs(g01_eqn, s, s_val);
    g01_eqn = subs(g01_eqn, h1, h1_val);
    g01_eqn = subs(g01_eqn, h2, h2_val);
    g01_eqn = subs(g01_eqn, h3, h3_val);
    g01_eqn = subs(g01_eqn, beta, beta_val);
    g01_eqn = subs(g01_eqn, gamma, gamma_val);

    g11_eqn = subs(mut_g11_eqn, mu, mu_val);
    g11_eqn = subs(g11_eqn, nu, nu_val);
    g11_eqn = subs(g11_eqn, s, s_val);
    g11_eqn = subs(g11_eqn, h1, h1_val);
    g11_eqn = subs(g11_eqn, h2, h2_val);
    g11_eqn = subs(g11_eqn, h3, h3_val);
    g11_eqn = subs(g11_eqn, beta, beta_val);
    g11_eqn = subs(g11_eqn, gamma, gamma_val);

    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end