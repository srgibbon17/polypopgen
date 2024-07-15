fig_HE2 = openfig("lower h2 2HEs extended.fig", 'invisible');

axes_HE2 = fig_HE2.Children;
data_HE2 = axes_HE2.Children;

x_HE2 = data_HE2.XData;
y_HE2 = data_HE2.YData;

fig_HE1 = openfig('lower h2 1HE extended.fig','invisible');

axes_HE1 = fig_HE1.Children;
data_HE1 = axes_HE1.Children;

x_HE1 = data_HE1.XData;
y_HE1 = data_HE1.YData;

fig_HE0 = openfig("lower h2 0HEs extended.fig", 'invisible');

axes_low_HE0 = fig_HE0.Children;
data_low_HE0 = axes_low_HE0.Children;

x_HE0 = data_low_HE0.XData;
y_HE0 = data_low_HE0.YData;

delta_HE2_HE1 = y_HE2 - y_HE1;
delta_HE2_HE0 = y_HE2 - y_HE0;
delta_HE1_HE0 = y_HE1 - y_HE0;

figure

title('Change in q across dominance (fixed mu)')

subplot(1, 3, 1)
scatter(x_HE2, delta_HE2_HE1);
xscale log
yscale log
title('HE2 - HE1')
ylabel('Change in q')

subplot(1, 3, 2)
scatter(x_HE2, delta_HE1_HE0)
xscale log
yscale log
title('HE1 - HE0')
xlabel('s')

subplot(1, 3, 3)
scatter(x_HE2, delta_HE2_HE0)
xscale log
yscale log
title('HE2 - HE0')

%subplot(1, 4, 4)
%scatter(x_HE2, delta_delta)
%xscale log
%yscale log
%title('Difference of differences')