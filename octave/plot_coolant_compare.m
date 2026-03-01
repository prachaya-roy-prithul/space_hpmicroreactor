% Coolant comparison plots (withdrawn vs inserted)
%
% Inputs:
% - results/coolant_compare.csv
%
% Outputs:
% - results/coolant_keff_compare.png
% - results/coolant_leakage_compare.png
% - results/coolant_shutdown_margin.png

clear; clc;

% ----------------------------
% Input path
% ----------------------------

in_csv = 'results/coolant_compare.csv';

% ----------------------------
% Read CSV (robust, no readtable)
% Columns:
% case,Tfuel_K,reflector_thickness_cm,rod_state,coolant,keff,unc,leak_frac,leak_frac_unc
% ----------------------------

fid = fopen(in_csv, 'r');
if fid < 0
  error('Could not open %s', in_csv);
end

C = textscan(fid, ...
  '%s %f %f %s %s %f %f %f %f', ...
  'Delimiter', ',', ...
  'HeaderLines', 1, ...
  'CollectOutput', false);

fclose(fid);

case_name = C{1};
Tfuel     = C{2};
refl_cm   = C{3};
rod_state = C{4};
coolant   = C{5};
keff      = C{6};
keff_unc  = C{7};
leak      = C{8};
leak_unc  = C{9};

n = length(keff);
if n == 0
  error('No data rows read from %s (file may be empty or delimiter mismatch).', in_csv);
end

% ----------------------------
% Sanity checks
% ----------------------------

Tf_unique = unique(Tfuel);
refl_unique = unique(refl_cm);

if length(Tf_unique) ~= 1
  fprintf('WARNING: Multiple Tfuel values detected. Using the first: %.1f K\n', Tf_unique(1));
end
if length(refl_unique) ~= 1
  fprintf('WARNING: Multiple reflector thickness values detected. Using the first: %.3f cm\n', refl_unique(1));
end

Tf_plot = Tf_unique(1);
refl_plot = refl_unique(1);

% ----------------------------
% Coolant ordering
% ----------------------------

order = {'Na', 'K', 'NaK', 'Li'};

present = false(size(order));
for j = 1:length(order)
  present(j) = any(strcmp(coolant, order{j}));
end
coolants = order(present);
m = length(coolants);

if m == 0
  error('No recognized coolant labels found. Expected one of: Na, K, NaK, Li.');
end

% ----------------------------
% Build aligned arrays (withdrawn vs inserted)
% ----------------------------

k_w  = nan(m,1);  u_w  = nan(m,1);
k_i  = nan(m,1);  u_i  = nan(m,1);

l_w  = nan(m,1);  lu_w = nan(m,1);
l_i  = nan(m,1);  lu_i = nan(m,1);

for i = 1:n
  cj = find(strcmp(coolants, coolant{i}));
  if isempty(cj)
    continue;
  end

  if strcmp(rod_state{i}, 'withdrawn')
    k_w(cj)  = keff(i);
    u_w(cj)  = keff_unc(i);
    l_w(cj)  = leak(i);
    lu_w(cj) = leak_unc(i);
  elseif strcmp(rod_state{i}, 'inserted')
    k_i(cj)  = keff(i);
    u_i(cj)  = keff_unc(i);
    l_i(cj)  = leak(i);
    lu_i(cj) = leak_unc(i);
  end
end

% ----------------------------
% Compute shutdown margin in pcm: (rho_w - rho_i)*1e5
% rho = (k - 1)/k
% Propagation: sigma_rho = sigma_k / k^2
% sigma_SM = 1e5 * sqrt(sigma_rho_w^2 + sigma_rho_i^2)
% ----------------------------

rho_w = (k_w - 1) ./ k_w;
rho_i = (k_i - 1) ./ k_i;

sig_rho_w = u_w ./ (k_w.^2);
sig_rho_i = u_i ./ (k_i.^2);

sm_pcm = (rho_w - rho_i) * 1.0e5;
sm_pcm_unc = 1.0e5 * sqrt(sig_rho_w.^2 + sig_rho_i.^2);

% ----------------------------
% Plot styling helpers
% ----------------------------

x = (1:m)';
xlab = coolants;

out_dir = 'results';
if exist(out_dir, 'dir') ~= 7
  mkdir(out_dir);
end

grid_alpha = 0.15;

% ----------------------------
% 1) k_eff vs coolant (withdrawn vs inserted)
% ----------------------------

figure('color', 'w');

hold on;

h1 = errorbar(x - 0.06, k_w, u_w, 'o');
set(h1, 'linewidth', 1.5);
set(h1, 'markersize', 8);

h2 = errorbar(x + 0.06, k_i, u_i, 's');
set(h2, 'linewidth', 1.5);
set(h2, 'markersize', 8);

% k = 1 reference line
plot([0.5 (m + 0.5)], [1 1], '--k', 'linewidth', 1.2);

hold off;

set(gca, 'xtick', x);
set(gca, 'xticklabel', xlab);

xlabel('Coolant', 'fontsize', 14);
ylabel('Effective multiplication factor, k_{eff}', 'fontsize', 14);

title({sprintf('k_{eff} comparison across coolants (withdrawn vs inserted)'), ...
       sprintf('T_{fuel} = %.0f K, reflector thickness = %.1f cm', Tf_plot, refl_plot)}, ...
      'fontsize', 14);

legend({'Withdrawn (void rod channel)', 'Inserted (B4C rod)', 'k = 1'}, 'location', 'northwest');
grid on;
set(gca, 'gridalpha', grid_alpha);
box on;

print('-dpng', '-r300', fullfile(out_dir, 'coolant_keff_compare.png'));
fprintf('Saved: %s\n', fullfile(out_dir, 'coolant_keff_compare.png'));

% ----------------------------
% 2) Leakage fraction vs coolant (withdrawn vs inserted)
% ----------------------------

figure('color', 'w');

hold on;

h1 = errorbar(x - 0.06, l_w, lu_w, 'o');
set(h1, 'linewidth', 1.5);
set(h1, 'markersize', 8);

h2 = errorbar(x + 0.06, l_i, lu_i, 's');
set(h2, 'linewidth', 1.5);
set(h2, 'markersize', 8);

hold off;

set(gca, 'xtick', x);
set(gca, 'xticklabel', xlab);

xlabel('Coolant', 'fontsize', 14);
ylabel('Leakage fraction (tally-derived)', 'fontsize', 14);

title({sprintf('Leakage fraction across coolants (withdrawn vs inserted)'), ...
       sprintf('T_{fuel} = %.0f K, reflector thickness = %.1f cm', Tf_plot, refl_plot)}, ...
      'fontsize', 14);

legend({'Withdrawn', 'Inserted'}, 'location', 'northeast');
grid on;
set(gca, 'gridalpha', grid_alpha);
box on;

print('-dpng', '-r300', fullfile(out_dir, 'coolant_leakage_compare.png'));
fprintf('Saved: %s\n', fullfile(out_dir, 'coolant_leakage_compare.png'));

% ----------------------------
% 3) Shutdown margin vs coolant (rho_w - rho_i), with propagated uncertainty
% ----------------------------

figure('color', 'w');

h = errorbar(x, sm_pcm, sm_pcm_unc, 'o');
set(h, 'linewidth', 1.5);
set(h, 'markersize', 8);

set(gca, 'xtick', x);
set(gca, 'xticklabel', xlab);

xlabel('Coolant', 'fontsize', 14);
ylabel('Shutdown margin, \rho_{w} - \rho_{i} (pcm)', 'fontsize', 14);

title({sprintf('Shutdown margin across coolants'), ...
       sprintf('T_{fuel} = %.0f K, reflector thickness = %.1f cm', Tf_plot, refl_plot)}, ...
      'fontsize', 16);

grid on;
set(gca, 'gridalpha', grid_alpha);
box on;

print('-dpng', '-r300', fullfile(out_dir, 'coolant_shutdown_margin.png'));
fprintf('Saved: %s\n', fullfile(out_dir, 'coolant_shutdown_margin.png'));