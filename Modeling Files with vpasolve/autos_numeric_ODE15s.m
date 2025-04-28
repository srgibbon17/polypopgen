function solve_ode_system()
    % Define parameter values
    mu_val = 2e-8;
    nu_val = 1e-9;
    a_val = 0;
    h1_val = 1;
    h2_val = 1;
    h3_val = 1;
    s_val = 1e-6;
    
    % Initial conditions (must sum to <= 1, as g2 = 1 - g0 - g1)
    g0_init = 0.04;
    g1_init = 0.32;
    
    % Time span
    tspan = [0 10^11];
    
    % Solve the ODE system
    [t, y] = ode15s(@(t, y) odefunc(t, y, s_val, mu_val, nu_val, a_val, h1_val, h2_val, h3_val), tspan, [g0_init, g1_init]);
    
    % Plot results
    figure;
    plot(t, y(:,1), 'r', 'DisplayName', 'g0'); hold on;
    plot(t, y(:,2), 'b', 'DisplayName', 'g1');
    plot(t, 1 - y(:,1) - y(:,2), 'g', 'DisplayName', 'g2');
    xlabel('Time'); ylabel('Frequencies');
    xscale log
    legend; grid on;
    title('Evolution of g0, g1, and g2');
end

function dydt = odefunc(~, y, s, mu, nu, a, h1, h2, h3)
    % Extract variables
    g0 = y(1);
    g1 = y(2);
    g2 = 1 - g0 - g1; % Constrained by simplex
    
    % Relative fitness calculations
    wbar = 1 - s * (2*g0*g1*h1 + (2*g0*g2 + g1^2)*h2 + 2*g1*g2*h3 + g2^2);
    w0 = (1/wbar);
    w1 = ((1 - s*h1)/wbar);
    w2 = ((1 - s*h2)/wbar);
    w3 = ((1 - s*h3)/wbar);
    w4 = ((1 - s)/wbar);
    
    % Selection terms
    sel_g0 = g0^2 * w0 + (1/2 + a/4) * (2*g0*g1) * w1 + (1/6 + a/3) * (2*g0*g2 + g1^2) * w2 + (a/4) * (2*g1*g2) * w3;
    sel_g1 = (1/2 - a/2) * (2*g0*g1) * w1 + (2/3 - 2*a/3) * (2*g0*g2 + g1^2) * w2 + (1/2 - a/2) * (2*g1*g2) * w3;
    
    % Mutation terms
    mut_g0 = sel_g0 * ((1 - mu)^2) + sel_g1 * (1 - mu) * nu + g2 * (nu^2) - g0;
    mut_g1 = 2 * sel_g0 * (1 - mu) * mu + sel_g1 * (1 - mu) * (1 - nu) + 2 * g2 * (1 - nu) * nu - g1;
    
    % ODE system
    dydt = [mut_g0; mut_g1];
end

solve_ode_system