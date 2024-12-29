syms c b beta gamma beta_prime

eqn_1 = c == (1-gamma)*(beta)^2 + (beta)*gamma;
eqn_2 = b == 2*(1-gamma)*(beta)*(beta_prime);
eqn_3 = 1-b-c == (1-gamma)*(beta_prime)^2+gamma*(beta_prime);
eqn_4 = beta == 1-beta_prime;

solve([eqn_1, eqn_2, eqn_3, eqn_4], gamma, "ReturnConditions", true)

eqn_1 = c == (1-gamma)*(b+c)^2 + (b+c)*gamma;
eqn_2 = b == 2*(1-gamma)*(b+c)*(1-b-c);
eqn_3 = 1-b-c == (1-gamma)*(1-b-c)^2+gamma*(1-b-c);

solve([eqn_1, eqn_2, eqn_3], gamma, "ReturnConditions", true)