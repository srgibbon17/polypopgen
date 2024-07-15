fig_2HE = openfig("2HEs 30 iterations.fig", 'invisible');

axes_2HE = fig_2HE.Children;
data_2HE = axes_2HE.Children;

x_2HE = data_2HE.XData;
y_2HE = data_2HE.YData;
z_2HE = data_2HE.ZData;

fig_1HE = openfig('1HE 30 iterations.fig','invisible');

axes_1HE = fig_1HE.Children;
data_1HE= axes_1HE.Children;

x_1HE = data_1HE.XData;
y_1HE = data_1HE.YData;
z_1HE = data_1HE.ZData;

fig_0HE = openfig("0HEs 30 iterations.fig", 'invisible');

axes_0HE = fig_0HE.Children;
data_0HE = axes_0HE.Children;

x_0HE = data_0HE.XData;
y_0HE = data_0HE.YData;
z_0HE = data_0HE.ZData;

delta_2HEs_1HE = z_2HE - z_1HE;
delta_2HEs_0HEs = z_2HE - z_0HE;
delta_1HE_0HEs = z_1HE - z_0HE;

figure

scatter3(x_2HE, y_2HE, delta_2HEs_1HE);
title('Change in q across HEs - 2HEs-1HE')
xscale log
yscale log
zscale log
xlabel('s')
ylabel('mu')
zlabel('Change in q')

figure

scatter3(x_2HE,  y_2HE, delta_1HE_0HEs);
title('Change in q across HEs - 1HE-0HEs')
xscale log
yscale log
zscale log
xlabel('s')
ylabel('mu')
zlabel('Change in q')

figure

scatter3(x_2HE,  y_2HE, delta_2HEs_0HEs);
title('Change in q across HEs - 2HEs-0HEs')
xscale log
yscale log
zscale log
xlabel('s')
ylabel('mu')
zlabel('Change in q')


%subplot(1, 4, 4)
%scatter(x_HE2, delta_delta)
%xscale log
%yscale log
%title('Difference of differences')