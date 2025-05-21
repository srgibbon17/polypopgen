function solve_ode_system()
    % Define parameter values
    mu_val = 2e-8;
    nu_val = 1e-9;
    a_val = 0;
    h1_val = 1;
    h2_val = 1;
    h3_val = 1;
    s_val = 1e-6;
    
    % Initial conditions for log-transformed variables
    g0_init = 0.04;
    g1_init = 0.32;
    g2_init = 1 - g0_init - g1_init;
    L0_init = log(g0_init);
    L1_init = log(g1_init);
    L2_init = log(g2_init);
    
    % Time span
    tspan = [0 10^10];
    
    % Solve the ODE system
    [t, y] = ode15s(@(t, y) odefunc(t, y, s_val, mu_val, nu_val, a_val, h1_val, h2_val, h3_val), tspan, [L0_init, L1_init, L2_init]);
    
    % Convert back to original variables
    g0 = exp(y(:,1));
    g1 = exp(y(:,2));
    g2 = exp(y(:,3));
    
    % Plot results
    figure;
    plot(t, g0, 'r', 'DisplayName', 'g0'); hold on;
    plot(t, g1, 'b', 'DisplayName', 'g1');
    plot(t, g2, 'g', 'DisplayName', 'g2');
    xlabel('Time'); ylabel('Frequencies');
    set(gca, 'XScale', 'log');
    legend; grid on;
    title('Evolution of g0, g1, and g2 (Log-Transformed Model)');
end

function dydt = odefunc(~, y, s, mu, nu, a, h1, h2, h3)
    % Extract log-transformed variables
    L0 = y(1);
    L1 = y(2);
    L2 = y(3);
    
    % Convert to original variables
    g0 = exp(L0);
    g1 = exp(L1);
    g2 = exp(L2);
    
    % % Ensure positivity
    % if g0 <= 0 || g1 <= 0 || g2 <= 0
    %     dydt = [0; 0; 0];
    %     return;
    % end
    
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
    mut_g2 = -mut_g0 - mut_g1;
    
    % Log-transformed ODEs
    dL0dt = (1 / g0) * mut_g0;
    dL1dt = (1 / g1) * mut_g1;
    dL2dt = (1 / g2) * mut_g2;
    
    dydt = [dL0dt; dL1dt; dL2dt];
end

solve_ode_system
