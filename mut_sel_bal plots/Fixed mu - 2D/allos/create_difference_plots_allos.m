fig_2HE = openfig("HE2 fixed mu zoomed.fig", 'invisible');

axes_2HE = fig_2HE.Children;
data_2HE = axes_2HE.Children;

x_2HE = data_2HE.XData;
y_2HE = data_2HE.YData;

fig_1HE = openfig('HE1 fixed mu zoomed.fig','invisible');

axes_1HE = fig_1HE.Children;
data_1HE= axes_1HE.Children;

x_1HE = data_1HE.XData;
y_1HE = data_1HE.YData;

fig_0HE = openfig("HE0 fixed mu zoomed.fig", 'invisible');

axes_0HE = fig_0HE.Children;
data_0HE = axes_0HE.Children;

x_0HE = data_0HE.XData;
y_0HE = data_0HE.YData;

delta_2HEs_1HE = y_2HE - y_1HE;
delta_2HEs_0HEs = y_2HE - y_0HE;
delta_1HE_0HEs = y_1HE - y_0HE;

figure

title('Change in q across alpha (fixed mu)')

subplot(1, 3, 1)
scatter(x_2HE, delta_2HEs_1HE);
xscale log
yscale log
title('2HEs - 1HE')
ylabel('Change in q')

subplot(1, 3, 2)
scatter(x_2HE, delta_1HE_0HEs)
xscale log
yscale log
title('1HE - 0HEs')
xlabel('s')

subplot(1, 3, 3)
scatter(x_2HE, delta_2HEs_0HEs)
xscale log
yscale log
title('2HEs - 0HEs')

%subplot(1, 4, 4)
%scatter(x_HE2, delta_delta)
%xscale log
%yscale log
%title('Difference of differences')