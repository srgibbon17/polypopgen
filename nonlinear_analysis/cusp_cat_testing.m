a_val = 0;
mu_val = 2e-8;
nu_val = 1e-9;
s_val_range = logspace(-9, -4, 50);
k_val_range = linspace(2/3, 1, 20);

[neutral_q, selected_q, unstable_q, neutral_avg_fitness, selected_avg_fitness, unstable_avg_fitness, s_coord, k_coord] = auto_cusp_cat_data(s_val_range, mu_val, nu_val, k_val_range, a_val);


h2_coord = 1-k_coord;


% k parameter without avg fitness coloring
figure

surf(s_coord, k_coord, unstable_q)

hold on

surf(s_coord, k_coord, selected_q)
surf(s_coord, k_coord, neutral_q)

xscale log
xlabel('s')
ylabel('k')
zlabel('q')
shading interp

% h_2 without avg fitness coloring
figure

surf(s_coord, h2_coord, unstable_q)

hold on

surf(s_coord, h2_coord, selected_q)
surf(s_coord, h2_coord, neutral_q)

xscale log
xlabel('s')
ylabel('h_2')
zlabel('q')
shading interp




