% Reflector sweep plot: k_eff vs reflector thickness

data = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);

thickness = data(:,2);
k         = data(:,3);
unc       = data(:,4);

% Sort just in case
[thickness, idx] = sort(thickness);
k   = k(idx);
unc = unc(idx);

figure('color','w');

% Plot with error bars
h = errorbar(thickness, k, unc, 'o-');
set(h, 'linewidth', 1.5);
set(h, 'markersize', 7);

hold on;

% Horizontal k = 1 line
yline = 1.0;
plot([min(thickness) max(thickness)], [yline yline], '--k', 'linewidth', 1.2);

% Estimate critical thickness by linear interpolation
% Find first crossing
for i = 1:length(k)-1
    if k(i) < 1 && k(i+1) > 1
        tcrit = thickness(i) + ...
            (1 - k(i)) * (thickness(i+1) - thickness(i)) / (k(i+1) - k(i));
        break;
    end
end

% Plot critical thickness marker
plot(tcrit, 1, 'sr', 'markersize', 8, 'markerfacecolor', 'r');

% Annotation
text(tcrit, 1.002, sprintf('  Critical thickness ≈ %.2f cm', tcrit), ...
     'fontsize', 11);

% Tight axis limits
xlim([0 12]);
ylim([0.70 1.05]);

xlabel('Reflector thickness (cm)', 'fontsize', 12);
ylabel('Effective multiplication factor, k_{eff}', 'fontsize', 12);

title({'Effective multiplication factor vs reflector thickness', ...
       'T_{fuel} = 600 K, control rod withdrawn'}, ...
       'fontsize', 13);

grid on;
set(gca, 'gridalpha', 0.15);  % Light grid

box on;

% Export high resolution
print('-dpng', '-r300', 'results/reflector_keff_withdrawn.png');

fprintf('Saved reflector_keff figure to results/reflector_keff_withdrawn.png\n');