%diploid data
dip_rec = readmatrix('dip_rec.csv');
dip_add = readmatrix('dip_add.csv');
dip_dom = readmatrix('dip_dom.csv');

%auto data
auto_rec = readmatrix('auto_rec.csv');
auto_add = readmatrix('auto_add.csv');
auto_dom = readmatrix('auto_dom.csv');

%allo data
allo_rec = readmatrix('allo_rec.csv');
allo_add = readmatrix('allo_add.csv');
allo_dom = readmatrix('allo_dom.csv');

dip_color = '#F04D13';
auto_color = '#66BED6';
allo_color = '#7DAB5B';

mu_val = 1e-8;
nu_val = 1e-8;
sel_range = dip_rec(:, 1);

%%%PLOTTTING%%%

figure

subplot(2, 3, 1)
% recessive allele frequency
plot(dip_rec(:, 1), dip_rec(:, 2), 'Color', dip_color, 'LineWidth', 3, 'DisplayName', 'Diploids')
hold on
plot(auto_rec(:, 1), auto_rec(:, 2), 'Color', auto_color, 'LineWidth', 3, 'DisplayName', ['Auto- and ' ...
    'Allotetraploids'])
%plot(allo_rec(:, 1), allo_rec(:, 2), 'Color', allo_color, 'LineWidth', 3, 'DisplayName', 'Allotetraploids')
plot(sel_range, (mu_val./sel_range).^(1/4), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2)
plot(sel_range, ((mu_val./sel_range) + (nu_val./(1-sel_range))).^(1/2), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 2)

xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Fully Recessive')

ylabel('q (Allele Frequency)')

subplot(2, 3, 2)
% additive allele frequency
plot(dip_add(:, 1), dip_add(:, 2), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(auto_add(:, 1), auto_add(:, 2), 'Color', auto_color, 'LineWidth', 3)
%plot(allo_add(:, 1), allo_add(:, 2), 'Color', allo_color, 'LineWidth', 3)

plot(sel_range, (mu_val./(sel_range.*.5)), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)
plot(sel_range, (mu_val./(sel_range.*.25)), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2)


xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Additive')
xlabel('s (Selection Coefficient)')

subplot(2, 3, 3)
% dominant allele frequency
plot(dip_dom(:, 1), dip_dom(:, 2), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(auto_dom(:, 1), auto_dom(:, 2), 'Color', auto_color, 'LineWidth', 3)
%plot(allo_dom(:, 1), allo_dom(:, 2), 'Color', allo_color, 'LineWidth', 3)

plot(sel_range, (mu_val./(sel_range)), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)
plot(sel_range, (mu_val./(sel_range)), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2)

xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Fully Dominant')

subplot(2, 3, 4)
% recessive average fitness
plot(dip_rec(:, 1), 1-dip_rec(:, 5), 'Color', dip_color, 'LineWidth', 3, 'DisplayName', 'Diploids')
hold on
plot(auto_rec(:, 1), 1-auto_rec(:, 6), 'Color', auto_color, 'LineWidth', 3, 'DisplayName', 'Autotetraploids')
%plot(allo_rec(:, 1), 1-allo_rec(:, 7), 'Color', allo_color, 'LineWidth', 3, 'DisplayName', 'Allotetraploids')
xscale log
xlim([1e-9, 1e-3])
ylim([0, 2e-8])

ylabel('Mutation Load')

subplot(2, 3, 5)
% additive average fitness
plot(dip_add(:, 1), 1-dip_add(:, 5), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(auto_add(:, 1), 1-auto_add(:, 6), 'Color', auto_color, 'LineWidth', 3)
%plot(allo_add(:, 1), 1-allo_add(:, 7), 'Color', allo_color, 'LineWidth', 3)

xscale log
xlim([1e-9, 1e-3])
xlabel('s (Selection Coefficient)')
ylim([0, 12e-8])



subplot(2, 3, 6)
% dominant average fitness
plot(dip_dom(:, 1), 1-dip_dom(:, 5), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(auto_dom(:, 1), 1-auto_dom(:, 6), 'Color', auto_color, 'LineWidth', 3)
%plot(allo_dom(:, 1), 1-allo_dom(:, 7), 'Color', allo_color, 'LineWidth', 3)
xscale log
xlim([1e-9, 1e-3])
ylim([0, 12e-8])





%%%%%%%%%%%%%%%
%Second Figure%
%   vertical  %
%%%%%%%%%%%%%%%

% figure
% 
% subplot(3, 2, 1)
% % recessive allele frequency
% plot(dip_rec(:, 1), dip_rec(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
% hold on
% plot(auto_rec(:, 1), auto_rec(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
% plot(allo_rec(:, 1), allo_rec(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')
% 
% xscale log
% xlim([1e-9, 1e-3])
% ylabel('Fully Recessive')
% 
% title('q (Allele Frequency)')
% legend
% 
% subplot(3, 2, 3)
% % additive allele frequency
% plot(dip_add(:, 1), dip_add(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
% hold on
% plot(auto_add(:, 1), auto_add(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
% plot(allo_add(:, 1), allo_add(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
% 
% xscale log
% xlim([1e-9, 1e-3])
% ylabel('Additive')
% 
% 
% subplot(3, 2, 5)
% % dominant allele frequency
% plot(dip_dom_neutral(:, 1), dip_dom_neutral(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
% hold on
% plot(dip_dom_selected(:, 1), dip_dom_selected(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
% plot(dip_dom_unstable(:, 1), dip_dom_unstable(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'LineStyle', '--')
% 
% plot(auto_dom_neutral(:, 1), auto_dom_neutral(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
% plot(auto_dom_selected(:, 1), auto_dom_selected(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
% plot(auto_dom_unstable(:, 1), auto_dom_unstable(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'LineStyle', '--')
% 
% plot(allo_dom_neutral(:, 1), allo_dom_neutral(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
% plot(allo_dom_selected(:, 1), allo_dom_selected(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
% plot(allo_dom_unstable(:, 1), allo_dom_unstable(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'LineStyle', '--')
% 
% xscale log
% xlim([1e-9, 1e-3])
% xlabel('s (Selection Coefficient)')
% ylabel('Fully Dominant')
% 
% subplot(3, 2, 2)
% % recessive average fitness
% plot(dip_rec(:, 1), -1+dip_rec(:, 5), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
% hold on
% plot(auto_rec(:, 1), -1+auto_rec(:, 6), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
% plot(allo_rec(:, 1), -1+allo_rec(:, 7), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')
% xscale log
% xlim([1e-9, 1e-3])
% %yscale log
% 
% title('Mutation Load')
% 
% subplot(3, 2, 4)
% % additive average fitness
% plot(dip_add(:, 1), -1+dip_add(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
% hold on
% plot(auto_add(:, 1), -1+auto_add(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
% plot(allo_add(:, 1), -1+allo_add(:, 7), 'Color', allo_color, 'LineWidth', 1.5)
% 
% xscale log
% xlim([1e-9, 1e-3])
% %yscale log
% 
% subplot(3, 2, 6)
% % dominant average fitness
% 
% %plot(auto_dom_neutral(:, 1), -1+auto_dom_neutral(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
% %hold on
% plot(auto_dom_selected(:, 1), -1+auto_dom_selected(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
% hold on
% 
% %plot(allo_dom_neutral(:, 1), -1+allo_dom_neutral(:, 7), 'Color', allo_color, 'LineWidth', 1.5)
% plot(allo_dom_selected(:, 1), -1+allo_dom_selected(:, 7), 'Color', allo_color, 'LineWidth', 1.5)
% 
% %plot(dip_dom_neutral(:, 1), -1+dip_dom_neutral(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
% plot(dip_dom_selected(:, 1), -1+dip_dom_selected(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
% 
% xscale log
% xlim([1e-9, 1e-3])
% xlabel('s (Selection Coefficient)')
% %yscale log
