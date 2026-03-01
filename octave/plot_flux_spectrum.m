% Fuel-region neutron flux spectrum
% Track-length flux tally extracted from OpenMC statepoint.
% Spectrum plotted on log-log axes.

clear;
clc;

csv_path = "results/flux_spectrum.csv";

if exist(csv_path, "file") != 2
  error(["CSV not found: ", csv_path]);
endif

data = dlmread(csv_path, ",", 1, 0);

E_mid = data(:, 3);
phi = data(:, 4);

% Normalize spectrum for shape comparison
phi_norm = phi ./ max(phi);

figure;
loglog(E_mid, phi_norm, "LineWidth", 1.5);
grid on;

xlabel("Energy (eV)");
ylabel("Flux (track-length estimate, normalized)");
title("Fuel-Region Neutron Flux Spectrum");

set(gca, "FontSize", 12);

out_path = "results/flux_spectrum_loglog.png";
print(out_path, "-dpng", "-r300");

disp(["Wrote plot: ", out_path]);