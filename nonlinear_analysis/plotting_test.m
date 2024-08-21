% plotting test for alpha

x_vals = linspace(0, 10, 100);

y_vals_1 = linspace(0, 3, 100);

y_vals_2 = linspace(3.5, 6.5, 100);

y_vals_3 = linspace(7, 10, 100);

figure 
plot(x_vals, y_vals_1, 'Color',[0 0.4470 0.7410])
hold on
plot(x_vals, y_vals_2, 'Color',[0 0.4470 0.7410 .6])
plot(x_vals, y_vals_3, 'Color',[0 0.4470 0.7410 .2])
