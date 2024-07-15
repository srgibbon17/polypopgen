fig_1_6 = openfig("Autos a1_6 zoomed.fig", 'invisible');

axes_1_6 = fig_1_6.Children;
data_1_6 = axes_1_6.Children;

x_1_6 = data_1_6.XData;
y_1_6 = data_1_6.YData;

fig_1_12 = openfig("Autos a1_12 zoomed.fig", 'invisible');

axes_1_12 = fig_1_12.Children;
data_1_12= axes_1_12.Children;

x_1_12 = data_1_12.XData;
y_1_12 = data_1_12.YData;

fig_0 = openfig("Autos a0 zoomed.fig", 'invisible');

axes_0 = fig_0.Children;
data_0 = axes_0.Children;

x_0 = data_0.XData;
y_0 = data_0.YData;

delta_1_6_1_12 = y_1_6 - y_1_12;
delta_1_6_0 = y_1_6 - y_0;
delta_1_12_0 = y_1_12 - y_0;

figure


title('Change in q across alpha (fixed mu)')



subplot(1, 3, 1)
scatter(x_1_6, delta_1_6_1_12);
xscale log
yscale log
title('alpha 1/6 - 1/12')
ylabel('Change in q')

subplot(1, 3, 2)
scatter(x_1_6, delta_1_12_0)
xscale log
yscale log
title('alpha 1/12 - 0')
xlabel('s')

subplot(1, 3, 3)
scatter(x_1_6, delta_1_6_0)
xscale log
yscale log
title('alpha 1/6 - 0')

%subplot(1, 4, 4)
%scatter(x_HE2, delta_delta)
%xscale log
%yscale log
%title('Difference of differences')