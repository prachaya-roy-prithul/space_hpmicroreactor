% Shutdown margin vs reflector thickness (pcm)
% Definition: SM(t) = (rho_withdrawn(t) - rho_inserted(t)) * 1e5

w = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);
i = dlmread('results/reflector_sweep_inserted_clean.csv', ',', 0, 0);

tw = w(:,2); kw = w(:,3); uw = w(:,4);
ti = i(:,2); ki = i(:,3); ui = i(:,4);

[tw, idxw] = sort(tw); kw = kw(idxw); uw = uw(idxw);
[ti, idxi] = sort(ti); ki = ki(idxi); ui = ui(idxi);

if length(tw) ~= length(ti) || any(abs(tw - ti) > 1.0e-9)
  error('Thickness grids do not match between withdrawn and inserted clean CSVs.');
end

rho_w = (kw - 1.0) ./ kw;
rho_i = (ki - 1.0) ./ ki;

sm_pcm = 1.0e5 * (rho_w - rho_i);

drhodk_w = 1.0 ./ (kw .^ 2);
drhodk_i = 1.0 ./ (ki .^ 2);

rho_w_sd = abs(drhodk_w) .* uw;
rho_i_sd = abs(drhodk_i) .* ui;

sm_sd_pcm = 1.0e5 * sqrt(rho_w_sd.^2 + rho_i_sd.^2);

figure('color','w');

h = errorbar(tw, sm_pcm, sm_sd_pcm, 'o-');
set(h, 'linewidth', 1.5, 'markersize', 7);

xlim([0 12]);

ymin = min(sm_pcm - sm_sd_pcm) - 200;
ymax = max(sm_pcm + sm_sd_pcm) + 200;
ylim([ymin ymax]);

xlabel('Reflector thickness (cm)', 'fontsize', 12);
ylabel('Shutdown margin, \rho_{w} - \rho_{i} (pcm)', 'fontsize', 12);

title({'Shutdown margin vs reflector thickness', ...
       'T_{fuel} = 600 K'}, 'fontsize', 13);

grid on;
set(gca, 'gridalpha', 0.12);
box on;

print('-dpng', '-r300', 'results/shutdown_margin_vs_thickness.png');
fprintf('Saved results/shutdown_margin_vs_thickness.png\n');