% k_eff vs reflector thickness (withdrawn vs inserted), with critical thickness markers

w = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);
i = dlmread('results/reflector_sweep_inserted_clean.csv', ',', 0, 0);

tw = w(:,2); kw = w(:,3); sw = w(:,4);
ti = i(:,2); ki = i(:,3); si = i(:,4);

% Sort
[tw, idxw] = sort(tw); kw = kw(idxw); sw = sw(idxw);
[ti, idxi] = sort(ti); ki = ki(idxi); si = si(idxi);

figure('color','w');

% Error bars
hw = errorbar(tw, kw, sw, 'o-'); hold on;
hi = errorbar(ti, ki, si, 's-');

set(hw, 'linewidth', 1.8, 'markersize', 8);
set(hi, 'linewidth', 1.8, 'markersize', 8);

% k = 1 line
plot([min([tw;ti]) max([tw;ti])], [1 1], '--k', 'linewidth', 1.2);

% Critical thickness (linear interpolation) for withdrawn
tcrit_w = NaN;
for n = 1:length(kw)-1
    if (kw(n) < 1 && kw(n+1) > 1) || (kw(n) > 1 && kw(n+1) < 1)
        tcrit_w = tw(n) + (1 - kw(n)) * (tw(n+1) - tw(n)) / (kw(n+1) - kw(n));
        break;
    end
end

% Critical thickness (linear interpolation) for inserted
tcrit_i = NaN;
for n = 1:length(ki)-1
    if (ki(n) < 1 && ki(n+1) > 1) || (ki(n) > 1 && ki(n+1) < 1)
        tcrit_i = ti(n) + (1 - ki(n)) * (ti(n+1) - ti(n)) / (ki(n+1) - ki(n));
        break;
    end
end

% Markers + vertical guide lines
if ~isnan(tcrit_w)
    plot(tcrit_w, 1, 'sr', 'markersize', 10, 'markerfacecolor', 'r');
    plot([tcrit_w tcrit_w], [min([kw;ki]) 1], ':k', 'linewidth', 1.0);
    text(tcrit_w, 1.002, sprintf('  Withdrawn t_{crit} = %.2f cm', tcrit_w), 'fontsize', 11);
end

if ~isnan(tcrit_i)
    plot(tcrit_i, 1, 'dr', 'markersize', 10, 'markerfacecolor', 'r');
    plot([tcrit_i tcrit_i], [min([kw;ki]) 1], ':k', 'linewidth', 1.0);
    text(tcrit_i, 0.992, sprintf('  Inserted t_{crit} = %.2f cm', tcrit_i), 'fontsize', 11);
end

xlabel('Reflector thickness (cm)', 'fontsize', 13);
ylabel('Effective multiplication factor, k_{eff}', 'fontsize', 13);

title({'k_{eff} vs reflector thickness (withdrawn vs inserted)', ...
       'T_{fuel} = 600 K'}, 'fontsize', 14);

legend({'Withdrawn (void rod channel)', 'Inserted (B4C rod)', 'k = 1'}, 'location', 'northwest');

grid on;
set(gca, 'gridalpha', 0.12);
box on;

xlim([0 max([tw;ti])]);
ylim([min([kw;ki]) - 0.02, max([kw;ki]) + 0.02]);

print('-dpng', '-r300', 'results/reflector_keff_compare.png');
fprintf('Saved: results/reflector_keff_compare.png\n');