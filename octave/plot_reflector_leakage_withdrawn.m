% Reflector sweep: leakage fraction vs reflector thickness

data = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);

thickness = data(:,2);
leak      = data(:,5);
leak_unc  = data(:,6);

figure('color','w');

h = errorbar(thickness, leak, leak_unc, 'o-');
set(h, 'linewidth', 1.5);
set(h, 'markersize', 7);

xlim([0 12]);
ylim([0.52 0.70]);

xlabel('Reflector thickness (cm)', 'fontsize', 12);
ylabel('Leakage fraction', 'fontsize', 12);

title({'Neutron leakage fraction vs reflector thickness', ...
       '(Derived from surface current and absorption tallies)'}, ...
       'fontsize', 13);

grid on;
set(gca, 'gridalpha', 0.15);

box on;

print('-dpng', '-r300', 'results/reflector_leakage_withdrawn.png');

fprintf('Saved reflector_leakage_withdrawn figure to results/reflector_leakage_withdrawn.png\n');