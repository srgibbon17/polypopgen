%diploid data
dip_rec_unstable = readmatrix('dip_rec_unstable.csv');
dip_rec_selected = readmatrix('dip_rec_selected.csv');
dip_rec_neutral = readmatrix('dip_rec_neutral.csv');

dip_rec = vertcat(dip_rec_neutral, dip_rec_selected);

dip_add_unstable = readmatrix('dip_add_unstable.csv');
dip_add_selected = readmatrix('dip_add_selected.csv');
dip_add_neutral = readmatrix('dip_add_neutral.csv');

dip_add = vertcat(dip_add_neutral, dip_add_selected);

dip_dom_unstable = readmatrix('dip_dom_unstable.csv');
dip_dom_selected = readmatrix('dip_dom_selected.csv');
dip_dom_neutral = readmatrix('dip_dom_neutral.csv');

%auto data
auto_rec_unstable = readmatrix('auto_rec_unstable.csv');
auto_rec_selected = readmatrix('auto_rec_selected.csv');
auto_rec_neutral = readmatrix('auto_rec_neutral.csv');

auto_rec = vertcat(auto_rec_neutral, auto_rec_selected);

auto_add_unstable = readmatrix('auto_add_unstable.csv');
auto_add_selected = readmatrix('auto_add_selected.csv');
auto_add_neutral = readmatrix('auto_add_neutral.csv');

auto_add = vertcat(auto_add_neutral, auto_add_selected);

auto_dom_unstable = readmatrix('auto_dom_unstable.csv');
auto_dom_selected = readmatrix('auto_dom_selected.csv');
auto_dom_neutral = readmatrix('auto_dom_neutral.csv');

%allo data
allo_rec_unstable = readmatrix('allo_rec_unstable.csv');
allo_rec_selected = readmatrix('allo_rec_selected.csv');
allo_rec_neutral = readmatrix('allo_rec_neutral.csv');

allo_rec = vertcat(allo_rec_neutral, allo_rec_selected);

allo_add_unstable = readmatrix('allo_add_unstable.csv');
allo_add_selected = readmatrix('allo_add_selected.csv');
allo_add_neutral = readmatrix('allo_add_neutral.csv');

allo_add = vertcat(allo_add_neutral, allo_add_selected);

allo_dom_unstable = readmatrix('allo_dom_unstable.csv');
allo_dom_selected = readmatrix('allo_dom_selected.csv');
allo_dom_neutral = readmatrix('allo_dom_neutral.csv');

dip_color = '#F04D13';
auto_color = '#66BED6';
allo_color = '#7DAB5B';

%%%PLOTTTING%%%

figure

subplot(1, 3, 1)
% recessive allele frequency
plot(dip_rec(:, 1), dip_rec(:, 2), 'Color', dip_color, 'LineWidth', 3, 'DisplayName', 'Diploids')
hold on
plot(auto_rec(:, 1), auto_rec(:, 2), 'Color', auto_color, 'LineWidth', 3, 'DisplayName', ['Auto- and ' ...
    'Allotetraploids'])
%plot(allo_rec(:, 1), allo_rec(:, 2), 'Color', allo_color, 'LineWidth', 3, 'DisplayName', 'Allotetraploids')

xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Fully Recessive')

ylabel('q (Allele Frequency)')
legend

subplot(1, 3, 2)
% additive allele frequency
plot(dip_add(:, 1), dip_add(:, 2), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(auto_add(:, 1), auto_add(:, 2), 'Color', auto_color, 'LineWidth', 3)
%plot(allo_add(:, 1), allo_add(:, 2), 'Color', allo_color, 'LineWidth', 3)

xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Additive')
xlabel('s (Selection Coefficient)')

subplot(1, 3, 3)
% dominant allele frequency
plot(dip_dom_neutral(:, 1), dip_dom_neutral(:, 2), 'Color', dip_color, 'LineWidth', 3)
hold on
plot(dip_dom_selected(:, 1), dip_dom_selected(:, 2), 'Color', dip_color, 'LineWidth', 3)
plot(dip_dom_unstable(:, 1), dip_dom_unstable(:, 2), 'Color', dip_color, 'LineWidth', 3, 'LineStyle', '--')

plot(auto_dom_neutral(:, 1), auto_dom_neutral(:, 2), 'Color', auto_color, 'LineWidth', 3)
plot(auto_dom_selected(:, 1), auto_dom_selected(:, 2), 'Color', auto_color, 'LineWidth', 3)
plot(auto_dom_unstable(:, 1), auto_dom_unstable(:, 2), 'Color', auto_color, 'LineWidth', 3, 'LineStyle', '--')

%plot(allo_dom_neutral(:, 1), allo_dom_neutral(:, 2), 'Color', allo_color, 'LineWidth', 3)
%plot(allo_dom_selected(:, 1), allo_dom_selected(:, 2), 'Color', allo_color, 'LineWidth', 3)
%plot(allo_dom_unstable(:, 1), allo_dom_unstable(:, 2), 'Color', allo_color, 'LineWidth', 3, 'LineStyle', '--')

xscale log
xlim([1e-9, 1e-3])
ylim([0, 1])

title('Fully Dominant')

% subplot(2, 3, 4)
% % recessive average fitness
% plot(dip_rec(:, 1), 1-dip_rec(:, 5), 'Color', dip_color, 'LineWidth', 3, 'DisplayName', 'Diploids')
% hold on
% plot(auto_rec(:, 1), 1-auto_rec(:, 6), 'Color', auto_color, 'LineWidth', 3, 'DisplayName', 'Autotetraploids')
% %plot(allo_rec(:, 1), 1-allo_rec(:, 7), 'Color', allo_color, 'LineWidth', 3, 'DisplayName', 'Allotetraploids')
% xscale log
% xlim([1e-9, 1e-3])
% ylim([0, 2e-8])



% ylabel('Mutation Load')

% subplot(2, 3, 5)
% % additive average fitness
% plot(dip_add(:, 1), 1-dip_add(:, 5), 'Color', dip_color, 'LineWidth', 3)
% hold on
% plot(auto_add(:, 1), 1-auto_add(:, 6), 'Color', auto_color, 'LineWidth', 3)
% %plot(allo_add(:, 1), 1-allo_add(:, 7), 'Color', allo_color, 'LineWidth', 3)
% 
% xscale log
% xlim([1e-9, 1e-3])
% xlabel('s (Selection Coefficient)')
% ylim([0, 12e-8])
% 
% 
% 
% subplot(2, 3, 6)
% % dominant average fitness
% 
% %plot(auto_dom_neutral(:, 1), 1-auto_dom_neutral(:, 6), 'Color', auto_color, 'LineWidth', 3)
% %hold on
% plot(auto_dom_selected(:, 1), 1-auto_dom_selected(:, 6), 'Color', auto_color, 'LineWidth', 3)
% hold on
% 
% %plot(allo_dom_neutral(:, 1), 1-allo_dom_neutral(:, 7), 'Color', allo_color, 'LineWidth', 3)
% %plot(allo_dom_selected(:, 1), 1-allo_dom_selected(:, 7), 'Color', allo_color, 'LineWidth', 3)
% 
% %plot(dip_dom_neutral(:, 1), 1-dip_dom_neutral(:, 5), 'Color', dip_color, 'LineWidth', 3)
% plot(dip_dom_selected(:, 1), 1-dip_dom_selected(:, 5), 'Color', dip_color, 'LineWidth', 3)
% 
% xscale log
% xlim([1e-9, 1e-3])
% ylim([0, 12e-8])





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
