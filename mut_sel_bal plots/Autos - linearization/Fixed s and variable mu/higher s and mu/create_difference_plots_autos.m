fig_8_48 = openfig("a8_48.fig", 'invisible');
axes_8_48 = fig_8_48.Children;
data_8_48 = axes_8_48.Children;
x_8_48 = data_8_48.XData;
y_8_48 = data_8_48.YData;

fig_7_48 = openfig("a7_48.fig", 'invisible');
axes_7_48 = fig_7_48.Children;
data_7_48 = axes_7_48.Children;
x_7_48 = data_7_48.XData;
y_7_48 = data_7_48.YData;

fig_6_48 = openfig("a6_48.fig", 'invisible');
axes_6_48 = fig_6_48.Children;
data_6_48 = axes_6_48.Children;
x_6_48 = data_6_48.XData;
y_6_48 = data_6_48.YData;

fig_5_48 = openfig("a5_48.fig", 'invisible');
axes_5_48 = fig_5_48.Children;
data_5_48 = axes_5_48.Children;
x_5_48 = data_5_48.XData;
y_5_48 = data_5_48.YData;

fig_4_48 = openfig("a4_48.fig", 'invisible');
axes_4_48 = fig_4_48.Children;
data_4_48 = axes_4_48.Children;
x_4_48 = data_4_48.XData;
y_4_48 = data_4_48.YData;

fig_3_48 = openfig("a3_48.fig", 'invisible');
axes_3_48 = fig_3_48.Children;
data_3_48 = axes_3_48.Children;
x_3_48 = data_3_48.XData;
y_3_48 = data_3_48.YData;

fig_2_48 = openfig("a2_48.fig", 'invisible');
axes_2_48 = fig_2_48.Children;
data_2_48 = axes_2_48.Children;
x_2_48 = data_2_48.XData;
y_2_48 = data_2_48.YData;

fig_1_48 = openfig("a1_48.fig", 'invisible');
axes_1_48 = fig_1_48.Children;
data_1_48 = axes_1_48.Children;
x_1_48 = data_1_48.XData;
y_1_48 = data_1_48.YData;


fig_0 = openfig("a0.fig", 'invisible');
axes_0 = fig_0.Children;
data_0 = axes_0.Children;
x_0 = data_0.XData;
y_0 = data_0.YData;

delta_1 = y_1_48 - y_0;
delta_2 = y_2_48 - y_1_48;
delta_3 = y_3_48 - y_2_48;
delta_4 = y_4_48 - y_3_48;
delta_5 = y_5_48 - y_4_48;
delta_6 = y_6_48 - y_5_48;
delta_7 = y_7_48 - y_6_48;
delta_8 = y_8_48 - y_7_48;

data_1 = fitlm(log10(x_8_48), log10(delta_1));
betahat_1 = data_1.Coefficients.Estimate(1);

data_2 = fitlm(log10(x_8_48), log10(delta_2));
betahat_2 = data_2.Coefficients.Estimate(1);

data_3 = fitlm(log10(x_8_48), log10(delta_3));
betahat_3 = data_3.Coefficients.Estimate(1);

data_4 = fitlm(log10(x_8_48), log10(delta_4));
betahat_4 = data_4.Coefficients.Estimate(1);

data_5 = fitlm(log10(x_8_48), log10(delta_5));
betahat_5 = data_5.Coefficients.Estimate(1);

data_6 = fitlm(log10(x_8_48), log10(delta_6));
betahat_6 = data_6.Coefficients.Estimate(1);

data_7 = fitlm(log10(x_8_48), log10(delta_7));
betahat_7 = data_7.Coefficients.Estimate(1);

data_8 = fitlm(log10(x_8_48), log10(delta_8));
betahat_8 = data_8.Coefficients.Estimate(1);

intercepts = [betahat_1, betahat_2, betahat_3, betahat_4, betahat_5, betahat_6, betahat_7, betahat_8];
alpha_values = 1/48:1/48:8/48;

fitlm(alpha_values, intercepts)

figure
scatter(alpha_values, intercepts)
title('log regression intercepts vs. alpha values')
xlabel('alpha values')
ylabel('intercept values')

figure

title('Change in q across alpha values (fixed mu)')

subplot(2, 4, 1)
scatter(x_8_48, delta_1);
xscale log
yscale log
ylabel('Change in q')

subplot(2, 4, 2)
scatter(x_8_48, delta_2);
xscale log
yscale log

subplot(2, 4, 3)
scatter(x_8_48, delta_3);
xscale log
yscale log

subplot(2, 4, 4)
scatter(x_8_48, delta_4);
xscale log
yscale log

subplot(2, 4, 5)
scatter(x_8_48, delta_5);
xscale log
yscale log
xlabel('mu')

subplot(2, 4, 6)
scatter(x_8_48, delta_6);
xscale log
yscale log

subplot(2, 4, 7)
scatter(x_8_48, delta_7);
xscale log
yscale log

subplot(2, 4, 8)
scatter(x_8_48, delta_8);
xscale log
yscale log
