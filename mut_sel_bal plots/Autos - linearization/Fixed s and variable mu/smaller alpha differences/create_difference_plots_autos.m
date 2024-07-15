fig_10_60 = openfig("a10_60.fig", 'invisible');
axes_10_60 = fig_10_60.Children;
data_10_60 = axes_10_60.Children;
x_10_60 = data_10_60.XData;
y_10_60 = data_10_60.YData;

fig_9_60 = openfig("a9_60.fig", 'invisible');
axes_9_60 = fig_9_60.Children;
data_9_60 = axes_9_60.Children;
x_9_60 = data_9_60.XData;
y_9_60 = data_9_60.YData;

fig_8_60 = openfig("a8_60.fig", 'invisible');
axes_8_60 = fig_8_60.Children;
data_8_60 = axes_8_60.Children;
x_8_60 = data_8_60.XData;
y_8_60 = data_8_60.YData;

fig_7_60 = openfig("a7_60.fig", 'invisible');
axes_7_60 = fig_7_60.Children;
data_7_60 = axes_7_60.Children;
x_7_60 = data_7_60.XData;
y_7_60 = data_7_60.YData;

fig_6_60 = openfig("a6_60.fig", 'invisible');
axes_6_60 = fig_6_60.Children;
data_6_60 = axes_6_60.Children;
x_6_60 = data_6_60.XData;
y_6_60 = data_6_60.YData;

fig_5_60 = openfig("a5_60.fig", 'invisible');
axes_5_60 = fig_5_60.Children;
data_5_60 = axes_5_60.Children;
x_5_60 = data_5_60.XData;
y_5_60 = data_5_60.YData;

fig_4_60 = openfig("a4_60.fig", 'invisible');
axes_4_60 = fig_4_60.Children;
data_4_60 = axes_4_60.Children;
x_4_60 = data_4_60.XData;
y_4_60 = data_4_60.YData;

fig_3_60 = openfig("a3_60.fig", 'invisible');
axes_3_60 = fig_3_60.Children;
data_3_60 = axes_3_60.Children;
x_3_60 = data_3_60.XData;
y_3_60 = data_3_60.YData;

fig_2_60 = openfig("a2_60.fig", 'invisible');
axes_2_60 = fig_2_60.Children;
data_2_60 = axes_2_60.Children;
x_2_60 = data_2_60.XData;
y_2_60 = data_2_60.YData;

fig_1_60 = openfig("a1_60.fig", 'invisible');
axes_1_60 = fig_1_60.Children;
data_1_60 = axes_1_60.Children;
x_1_60 = data_1_60.XData;
y_1_60 = data_1_60.YData;


fig_0 = openfig("a0.fig", 'invisible');
axes_0 = fig_0.Children;
data_0 = axes_0.Children;
x_0 = data_0.XData;
y_0 = data_0.YData;

delta_1 = y_1_60 - y_0;
delta_2 = y_2_60 - y_1_60;
delta_3 = y_3_60 - y_2_60;
delta_4 = y_4_60 - y_3_60;
delta_5 = y_5_60 - y_4_60;
delta_6 = y_6_60 - y_5_60;
delta_7 = y_7_60 - y_6_60;
delta_8 = y_8_60 - y_7_60;
delta_9 = y_9_60 - y_8_60;
delta_10 = y_10_60 - y_9_60;

data_1 = fitlm(log10(x_8_60), log10(delta_1));
betahat_1 = data_1.Coefficients.Estimate(1);

data_2 = fitlm(log10(x_8_60), log10(delta_2));
betahat_2 = data_2.Coefficients.Estimate(1);

data_3 = fitlm(log10(x_8_60), log10(delta_3));
betahat_3 = data_3.Coefficients.Estimate(1);

data_4 = fitlm(log10(x_8_60), log10(delta_4));
betahat_4 = data_4.Coefficients.Estimate(1);

data_5 = fitlm(log10(x_8_60), log10(delta_5));
betahat_5 = data_5.Coefficients.Estimate(1);

data_6 = fitlm(log10(x_8_60), log10(delta_6));
betahat_6 = data_6.Coefficients.Estimate(1);

data_7 = fitlm(log10(x_8_60), log10(delta_7));
betahat_7 = data_7.Coefficients.Estimate(1);

data_8 = fitlm(log10(x_8_60), log10(delta_8));
betahat_8 = data_8.Coefficients.Estimate(1);

data_9 = fitlm(log10(x_8_60), log10(delta_9));
betahat_9 = data_9.Coefficients.Estimate(1);

data_10 = fitlm(log10(x_8_60), log10(delta_10));
betahat_10 = data_10.Coefficients.Estimate(1);

intercepts = [betahat_1, betahat_2, betahat_3, betahat_4, betahat_5, betahat_6, betahat_7, betahat_8, betahat_9, betahat_10];
alpha_values = 1/60:1/60:10/60;

fitlm(alpha_values, intercepts)

figure
scatter(alpha_values, intercepts)
title('log regression intercepts vs. alpha values')
xlabel('alpha values')
ylabel('intercept values')

figure

title('Change in q across alpha values (fixed mu)')

subplot(2, 5, 1)
scatter(x_8_60, delta_1);
xscale log
yscale log
ylabel('Change in q')

subplot(2, 5, 2)
scatter(x_8_60, delta_2);
xscale log
yscale log

subplot(2, 5, 3)
scatter(x_8_60, delta_3);
xscale log
yscale log

subplot(2, 5, 4)
scatter(x_8_60, delta_4);
xscale log
yscale log

subplot(2, 5, 5)
scatter(x_8_60, delta_5);
xscale log
yscale log
xlabel('mu')

subplot(2, 5, 6)
scatter(x_8_60, delta_6);
xscale log
yscale log

subplot(2, 5, 7)
scatter(x_8_60, delta_7);
xscale log
yscale log

subplot(2, 5, 8)
scatter(x_8_60, delta_8);
xscale log
yscale log

subplot(2, 5, 9)
scatter(x_8_60, delta_9);
xscale log
yscale log

subplot(2, 5, 10)
scatter(x_8_60, delta_10);
xscale log
yscale log
