
a_val = 0;
mu_val = 2e-8;
nu_val = 1e-9;
%s_val_range = logspace(-9, -4, 50);
%k_val_range = linspace(0, .1, 10);

% alternative s and k ranges with variably dense sampling:
s_val_range = [logspace(-9, -8, 11), logspace(-8, -5, 91), logspace(-5,-4, 11)];
k_val_range = linspace(0, 1, 51);

[neutral_q, selected_q, unstable_q, neutral_avg_fitness, selected_avg_fitness, unstable_avg_fitness, s_coord, k_coord] = auto_cusp_cat_data_TEST(s_val_range, mu_val, nu_val, k_val_range, a_val);


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

savefig('auto_cuspcat_k_dense_s.fig')

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

savefig('auto_cuspcat_h2_dense_s.fig')

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

savefig('auto_cuspcat_k_no_shading_dense_s.fig')

% plot with avg. fitness as colorbar
% figure
% 
% surf(s_coord, h2_coord, unstable_q, unstable_avg_fitness)
% 
% hold on
% 
% surf(s_coord, h2_coord, selected_q, selected_avg_fitness)
% surf(s_coord, h2_coord, neutral_q, neutral_avg_fitness)
% 
% xscale log
% xlabel('s')
% ylabel('h_2')
% zlabel('q')

%savefig('auto_cuspcat_h2_no_shading.fig')