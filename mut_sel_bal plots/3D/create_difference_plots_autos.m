fig_1_6 = openfig("autos a_1_6 recessive.fig", 'invisible');

axes_1_6 = fig_1_6.Children;
data_1_6 = axes_1_6.Children;

x_1_6 = data_1_6.XData;
y_1_6 = data_1_6.YData;
z_1_6 = data_1_6.ZData;

fig_1_12 = openfig('autos a_1_12 recessive.fig','invisible');

axes_1_12 = fig_1_12.Children;
data_1_12= axes_1_12.Children;

x_1_12 = data_1_12.XData;
y_1_12 = data_1_12.YData;
z_1_12 = data_1_12.ZData;

fig_0 = openfig("autos a_0 recessive.fig", 'invisible');

axes_0 = fig_0.Children;
data_0 = axes_0.Children;

x_0 = data_0.XData;
y_0 = data_0.YData;
z_0 = data_0.ZData;

delta_1_6_1_12 = z_1_6 - z_1_12;
delta_1_6_0 = z_1_6 - z_0;
delta_1_12_0 = z_1_12 - z_0;

figure

scatter3(x_1_6, y_1_6, delta_1_6_1_12);
title('Change in q across alpha - 1/6 - 1/12')
xscale log
yscale log
zscale log
xlabel('s')
ylabel('mu')
zlabel('Change in q')

figure

scatter3(x_1_6,  y_1_6, delta_1_12_0);
title('Change in q across HEs - 1/12 - 0')
xscale log
yscale log
zscale log
xlabel('s')
ylabel('mu')
zlabel('Change in q')

figure

scatter3(x_1_6,  y_1_6, delta_1_6_0);
title('Change in q across HEs - 1/6 - 0')
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