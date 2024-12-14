%diploid data
dip_rec_unstable = readmatrix('dip_rec_unstable.csv');
dip_rec_selected = readmatrix('dip_rec_selected.csv');
dip_rec_neutral = readmatrix('dip_rec_neutral.csv');

dip_add_unstable = readmatrix('dip_add_unstable.csv');
dip_add_selected = readmatrix('dip_add_selected.csv');
dip_add_neutral = readmatrix('dip_add_neutral.csv');

dip_add = vertcat(dip_add_neutral, dip_add_selected);

dip_dom_unstable = readmatrix('dip_dom_unstable.csv');
dip_dom_selected = readmatrix('dip_dom_selected.csv');
dip_dom_neutral = readmatrix('dip_dom_neutral.csv');

dip_dom = vertcat(dip_dom_neutral, dip_dom_selected);

%auto data
auto_rec_unstable = readmatrix('auto_rec_unstable.csv');
auto_rec_selected = readmatrix('auto_rec_selected.csv');
auto_rec_neutral = readmatrix('auto_rec_neutral.csv');

auto_add_unstable = readmatrix('auto_add_unstable.csv');
auto_add_selected = readmatrix('auto_add_selected.csv');
auto_add_neutral = readmatrix('auto_add_neutral.csv');

auto_add = vertcat(auto_add_neutral, auto_add_selected);

auto_dom_unstable = readmatrix('auto_dom_unstable.csv');
auto_dom_selected = readmatrix('auto_dom_selected.csv');
auto_dom_neutral = readmatrix('auto_dom_neutral.csv');

auto_dom = vertcat(auto_dom_neutral, auto_dom_selected);

%allo data
allo_rec_unstable = readmatrix('allo_rec_unstable.csv');
allo_rec_selected = readmatrix('allo_rec_selected.csv');
allo_rec_neutral = readmatrix('allo_rec_neutral.csv');

allo_add_unstable = readmatrix('allo_add_unstable.csv');
allo_add_selected = readmatrix('allo_add_selected.csv');
allo_add_neutral = readmatrix('allo_add_neutral.csv');

allo_add = vertcat(allo_add_neutral, allo_add_selected);

allo_dom_unstable = readmatrix('allo_dom_unstable.csv');
allo_dom_selected = readmatrix('allo_dom_selected.csv');
allo_dom_neutral = readmatrix('allo_dom_neutral.csv');

allo_dom = vertcat(allo_dom_neutral, allo_dom_selected);

dip_color = '#F04D13';
auto_color = '#66BED6';
allo_color = '#7DAB5B';

%%%PLOTTTING%%%

figure

subplot(2, 3, 1)
% recessive allele frequency
plot(dip_rec_neutral(:, 1), dip_rec_neutral(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(dip_rec_selected(:, 1), dip_rec_selected(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
plot(dip_rec_unstable(:, 1), dip_rec_unstable(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'LineStyle', '--')

plot(auto_rec_neutral(:, 1), auto_rec_neutral(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(auto_rec_selected(:, 1), auto_rec_selected(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(auto_rec_unstable(:, 1), auto_rec_unstable(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'LineStyle', '--')

plot(allo_rec_neutral(:, 1), allo_rec_neutral(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
plot(allo_rec_selected(:, 1), allo_rec_selected(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
plot(allo_rec_unstable(:, 1), allo_rec_unstable(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'LineStyle', '--')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

title('Fully Recessive')

ylabel('q (Derived Allele Frequency)')

subplot(2, 3, 2)
% additive allele frequency
plot(dip_add(:, 1), dip_add(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(auto_add(:, 1), auto_add(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(allo_add(:, 1), allo_add(:, 2), 'Color', allo_color, 'LineWidth', 1.5)

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

title('Additive')

subplot(2, 3, 3)
% dominant allele frequency
plot(dip_dom(:, 1), dip_dom(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
hold on
plot(auto_dom(:, 1), auto_dom(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
plot(allo_dom(:, 1), allo_dom(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

legend

title('Fully Dominant')

subplot(2, 3, 4)
% recessive average fitness

plot(auto_rec_neutral(:, 1), (auto_rec_neutral(:, 6)-1)./auto_rec_neutral(:,6), 'Color', auto_color, 'LineWidth', 1.5)
hold on
plot(auto_rec_selected(:, 1), (auto_rec_selected(:, 6)-1)./auto_rec_selected(:,6), 'Color', auto_color, 'LineWidth', 1.5)
%hold on

%plot(allo_rec_neutral(:, 1), (allo_rec_neutral(:, 7)-1)./allo_rec_neutral(:, 7), 'Color', allo_color, 'LineWidth', 1.5)
%plot(allo_rec_selected(:, 1), (allo_rec_selected(:, 7)-1)./allo_rec_selected(:, 7), 'Color', allo_color, 'LineWidth', 1.5)

plot(dip_rec_neutral(:, 1), (dip_rec_neutral(:, 5)-1)./dip_rec_neutral(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
plot(dip_rec_selected(:, 1), (dip_rec_selected(:, 5)-1)./dip_rec_selected(:, 5), 'Color', dip_color, 'LineWidth', 1.5)

plot(allo_rec_neutral(:, 1), -allo_rec_neutral(:,1)/(1+allo_rec_neutral(:,1)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')
plot(allo_rec_selected(:, 1), -allo_rec_selected(:,1)/(1+allo_rec_selected(:,1)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')
yscale log

ylabel('Deviation from Fitness Optimum')

subplot(2, 3, 5)
% additive average fitness
plot(dip_add(:, 1), (dip_add(:, 5)-1)./dip_add(:,5), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(auto_add(:, 1), (auto_add(:, 6)-1)./auto_add(:,6), 'Color', auto_color, 'LineWidth', 1.5)
plot(allo_add(:, 1), (allo_add(:, 7)-1)./allo_add(:,7), 'Color', allo_color, 'LineWidth', 1.5)

plot(allo_add(:, 1), -allo_add(:,1)/(1+allo_add(:,1)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')

xscale log
xlim([-1e-3, -1e-9])
xlabel('s (Selection Coefficient)')
set(gca, 'xdir', 'reverse')
yscale log

subplot(2, 3, 6)
% dominant average fitness

plot(dip_dom(:, 1), (dip_dom(:, 5)-1)./dip_dom(:,5), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
hold on
plot(auto_dom(:, 1), (auto_dom(:, 6)-1)./auto_dom(:,6), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
plot(allo_dom(:, 1), (allo_dom(:, 7)-1)./allo_dom(:,7), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')

plot(allo_dom(:, 1), -allo_dom(:,1)/(1+allo_dom(:,1)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')
yscale log

%%%%%%%%%%%%%%%
%Second Figure%
%%%%%%%%%%%%%%%

figure

subplot(3, 2, 1)
% recessive allele frequency
plot(dip_rec_neutral(:, 1), dip_rec_neutral(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(dip_rec_selected(:, 1), dip_rec_selected(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
plot(dip_rec_unstable(:, 1), dip_rec_unstable(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'LineStyle', '--')

plot(auto_rec_neutral(:, 1), auto_rec_neutral(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(auto_rec_selected(:, 1), auto_rec_selected(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(auto_rec_unstable(:, 1), auto_rec_unstable(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'LineStyle', '--')

plot(allo_rec_neutral(:, 1), allo_rec_neutral(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
plot(allo_rec_selected(:, 1), allo_rec_selected(:, 2), 'Color', allo_color, 'LineWidth', 1.5)
plot(allo_rec_unstable(:, 1), allo_rec_unstable(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'LineStyle', '--')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

title('q (Derived Allele Frequency)')
ylabel('Fully Recessive')

subplot(3, 2, 3)
% additive allele frequency
plot(dip_add(:, 1), dip_add(:, 2), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(auto_add(:, 1), auto_add(:, 2), 'Color', auto_color, 'LineWidth', 1.5)
plot(allo_add(:, 1), allo_add(:, 2), 'Color', allo_color, 'LineWidth', 1.5)

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

ylabel('Additive')

subplot(3, 2, 5)
% dominant allele frequency
plot(dip_dom(:, 1), dip_dom(:, 2), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
hold on
plot(auto_dom(:, 1), auto_dom(:, 2), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
plot(allo_dom(:, 1), allo_dom(:, 2), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')

legend

ylabel('Fully Dominant')

subplot(3, 2, 2)
% recessive average fitness

plot(auto_rec_neutral(:, 1), -1+auto_rec_neutral(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
hold on
plot(auto_rec_selected(:, 1), -1+auto_rec_selected(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
%hold on

plot(allo_rec_neutral(:, 1), -1+allo_rec_neutral(:, 7), 'Color', allo_color, 'LineWidth', 1.5)
plot(allo_rec_selected(:, 1), -1+allo_rec_selected(:, 7), 'Color', allo_color, 'LineWidth', 1.5)

plot(dip_rec_neutral(:, 1), -1+dip_rec_neutral(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
plot(dip_rec_selected(:, 1), -1+dip_rec_selected(:, 5), 'Color', dip_color, 'LineWidth', 1.5)

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')
yscale log

title('Deviation from Fitness Optimum')

subplot(3, 2, 4)
% additive average fitness
plot(dip_add(:, 1), -1+dip_add(:, 5), 'Color', dip_color, 'LineWidth', 1.5)
hold on
plot(auto_add(:, 1), -1+auto_add(:, 6), 'Color', auto_color, 'LineWidth', 1.5)
plot(allo_add(:, 1), -1+allo_add(:, 7), 'Color', allo_color, 'LineWidth', 1.5)

xscale log
xlim([-1e-3, -1e-9])
xlabel('s (Selection Coefficient)')
set(gca, 'xdir', 'reverse')
yscale log

subplot(3, 2, 6)
% dominant average fitness

plot(dip_dom(:, 1), -1+dip_dom(:, 5), 'Color', dip_color, 'LineWidth', 1.5, 'DisplayName', 'Diploids')
hold on
plot(auto_dom(:, 1), -1+auto_dom(:, 6), 'Color', auto_color, 'LineWidth', 1.5, 'DisplayName', 'Autotetraploids')
plot(allo_dom(:, 1), -1+allo_dom(:, 7), 'Color', allo_color, 'LineWidth', 1.5, 'DisplayName', 'Allotetraploids')

xscale log
xlim([-1e-3, -1e-9])
set(gca, 'xdir', 'reverse')
yscale log