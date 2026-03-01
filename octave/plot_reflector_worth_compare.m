% Reflector worth comparison plot (pcm): withdrawn vs inserted
% Worth definition: Delta rho(t) = (rho(t) - rho(0)) * 1e5

w = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);
i = dlmread('results/reflector_sweep_inserted_clean.csv', ',', 0, 0);

tw = w(:,2); kw = w(:,3); uw = w(:,4);
ti = i(:,2); ki = i(:,3); ui = i(:,4);

[tw, idxw] = sort(tw); kw = kw(idxw); uw = uw(idxw);
[ti, idxi] = sort(ti); ki = ki(idxi); ui = ui(idxi);

rho_w = (kw - 1.0) ./ kw;
rho_i = (ki - 1.0) ./ ki;

rho_w0 = rho_w(1);
rho_i0 = rho_i(1);

worth_w = 1.0e5 * (rho_w - rho_w0);
worth_i = 1.0e5 * (rho_i - rho_i0);

drhodk_w = 1.0 ./ (kw .^ 2);
drhodk_i = 1.0 ./ (ki .^ 2);

rho_w_sd = abs(drhodk_w) .* uw;
rho_i_sd = abs(drhodk_i) .* ui;

rho_w0_sd = rho_w_sd(1);
rho_i0_sd = rho_i_sd(1);

worth_w_sd = 1.0e5 * sqrt(rho_w_sd.^2 + rho_w0_sd.^2);
worth_i_sd = 1.0e5 * sqrt(rho_i_sd.^2 + rho_i0_sd.^2);

figure('color','w');
hold on;

hw = errorbar(tw, worth_w, worth_w_sd, 'o-');
set(hw, 'linewidth', 1.5, 'markersize', 7);

hi = errorbar(ti, worth_i, worth_i_sd, 's-');
set(hi, 'linewidth', 1.5, 'markersize', 7);

xlim([0 12]);

ymin = min([worth_w - worth_w_sd; worth_i - worth_i_sd]) - 500;
ymax = max([worth_w + worth_w_sd; worth_i + worth_i_sd]) + 500;
ylim([ymin ymax]);

xlabel('Reflector thickness (cm)', 'fontsize', 12);
ylabel('Reflector worth, \Delta\rho (pcm)', 'fontsize', 12);

title({'Reflector worth vs thickness (withdrawn vs inserted)', ...
       'Baseline: same rod state at 0 cm'}, 'fontsize', 13);

legend({'Withdrawn worth', 'Inserted worth'}, 'location', 'northwest');

grid on;
set(gca, 'gridalpha', 0.12);
box on;

print('-dpng', '-r300', 'results/reflector_worth_compare.png');
fprintf('Saved results/reflector_worth_compare.png\n');