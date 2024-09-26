
a_val = 1/6;
mu_val = 2e-8;
nu_val = 1e-9;

s_val_range = [logspace(-9, -8, 3), logspace(-8, -5.5, 14), logspace(-5.5,-4, 3)];
k_val_range = [linspace(0, .05, 20), linspace(.05, .2, 10), linspace(.2, 1, 5)];

k_bistable_range = [linspace(0, .05, 20), linspace(.05, .2, 10)];

[updated_s_val_range] = updated_s_grid(s_val_range, mu_val, nu_val, k_bistable_range, a_val);

disp('New s grid complete.')

[neutral_q, selected_q, unstable_q, neutral_avg_fitness, selected_avg_fitness, unstable_avg_fitness, s_coord, k_coord] = auto_cusp_cat_no_extrapolation(updated_s_val_range, mu_val, nu_val, k_val_range, a_val);


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

%savefig('auto_cuspcat_k.fig')

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
%shading interp

%savefig('auto_cuspcat_h2.fig')

% k param without shading
figure

surf(s_coord, k_coord, unstable_q)

hold on

surf(s_coord, k_coord, selected_q)
surf(s_coord, k_coord, neutral_q)

xscale log
xlabel('s')
ylabel('k')
zlabel('q')

%savefig('auto_cuspcat_k_no_shading.fig')

% plot with avg. fitness as colorbar
figure

surf(s_coord, h2_coord, unstable_q, unstable_avg_fitness)

hold on

surf(s_coord, h2_coord, selected_q, selected_avg_fitness)
surf(s_coord, h2_coord, neutral_q, neutral_avg_fitness)

xscale log
xlabel('s')
ylabel('h_2')
zlabel('q')

%savefig('auto_cuspcat_h2_no_shading.fig')