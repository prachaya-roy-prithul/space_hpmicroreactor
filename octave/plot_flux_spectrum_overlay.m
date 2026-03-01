% Fuel-region neutron flux spectrum overlay (withdrawn vs inserted).
% Each spectrum is normalized by its own maximum to compare spectral shape.

clear;
clc;

csv_w = "results/flux_spectrum_withdrawn.csv";
csv_i = "results/flux_spectrum_inserted.csv";

if exist(csv_w, "file") != 2
  error(["CSV not found: ", csv_w]);
endif

if exist(csv_i, "file") != 2
  error(["CSV not found: ", csv_i]);
endif

w = dlmread(csv_w, ",", 1, 0);
i = dlmread(csv_i, ",", 1, 0);

Ew = w(:, 3);
phiw = w(:, 4);

Ei = i(:, 3);
phii = i(:, 4);

phiw_n = phiw ./ max(phiw);
phii_n = phii ./ max(phii);

figure;
loglog(Ew, phiw_n, "LineWidth", 1.8);
hold on;
loglog(Ei, phii_n, "LineWidth", 1.8);
grid on;

xlabel("Energy (eV)");
ylabel("Flux (track-length estimate, normalized)");
title("Fuel-Region Neutron Flux Spectrum (Withdrawn vs Inserted)");

legend("Rod withdrawn", "Rod inserted", "location", "northwest");

set(gca, "FontSize", 12);

out_path = "results/flux_spectrum_overlay.png";
print(out_path, "-dpng", "-r300");
disp(["Wrote plot: ", out_path]);