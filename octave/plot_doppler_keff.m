% Doppler sweep 

% Outputs:
% - results/doppler_keff_errorbars.png
% - results/doppler_reactivity_errorbars.png
% - regression coefficients in pcm/K

clear;
clc;

csv_path = "results/doppler_keff.csv";

if exist(csv_path, "file") != 2
  error(["CSV not found: ", csv_path]);
endif

data = dlmread(csv_path, ",", 1, 0);

T = data(:, 2);
k = data(:, 3);
k_unc = data(:, 4);

[Ts, idx] = sort(T);
ks = k(idx);
k_unc_s = k_unc(idx);

% Reactivity and uncertainty propagation (first-order)
rho = (ks - 1.0) ./ ks;
drho_dk = 1.0 ./ (ks .^ 2);
rho_unc = abs(drho_dk) .* k_unc_s;

rho_pcm = rho * 1.0e5;
rho_unc_pcm = rho_unc * 1.0e5;

% Linear regression: k(T) and rho(T)
p_k = polyfit(Ts, ks, 1);
k_fit = polyval(p_k, Ts);
dkdT = p_k(1);
dkdT_pcm = dkdT * 1.0e5;

p_rho = polyfit(Ts, rho_pcm, 1);
rho_fit = polyval(p_rho, Ts);
drhodT_pcm = p_rho(1);

% Optional: sqrt(T) regression for Doppler-style scaling
x = sqrt(Ts);
p_rho_sqrt = polyfit(x, rho_pcm, 1);
rho_fit_sqrt = polyval(p_rho_sqrt, x);
drho_dsqrtT_pcm = p_rho_sqrt(1);

% ----------------------------
% Plot k vs T with error bars
% ----------------------------

figure;
h1 = errorbar(Ts, ks, k_unc_s);
set(h1, "marker", "o");
set(h1, "linestyle", "-");
set(h1, "linewidth", 1.3);

hold on;
h2 = plot(Ts, k_fit, "--");
set(h2, "linewidth", 1.3);

grid on;
xlabel("Fuel temperature (K)");
ylabel("k-effective");
title("Doppler Sweep: k-effective vs Fuel Temperature (with uncertainty)");
legend("OpenMC results", "Linear fit", "location", "southwest");
set(gca, "FontSize", 12);

print("results/doppler_keff_errorbars.png", "-dpng", "-r300");

% ----------------------------
% Plot reactivity vs T with error bars
% ----------------------------

figure;
h1 = errorbar(Ts, rho_pcm, rho_unc_pcm);
set(h1, "marker", "o");
set(h1, "linestyle", "-");
set(h1, "linewidth", 1.3);

hold on;
h2 = plot(Ts, rho_fit, "--");
set(h2, "linewidth", 1.3);

grid on;
xlabel("Fuel temperature (K)");
ylabel("Reactivity (pcm)");
title("Doppler Sweep: Reactivity vs Fuel Temperature (with uncertainty)");
legend("OpenMC results", "Linear fit", "location", "southwest");
set(gca, "FontSize", 12);

print("results/doppler_reactivity_errorbars.png", "-dpng", "-r300");

% ----------------------------
% Plot: Reactivity vs sqrt(T)
% ----------------------------

figure;
h1 = errorbar(x, rho_pcm, rho_unc_pcm);
set(h1, "marker", "o");
set(h1, "linestyle", "-");
set(h1, "linewidth", 1.3);

hold on;
h2 = plot(x, rho_fit_sqrt, "--");
set(h2, "linewidth", 1.3);

grid on;
xlabel("sqrt(Fuel temperature) (sqrt(K))");
ylabel("Reactivity (pcm)");
title("Doppler Sweep: Reactivity vs sqrt(Temperature)");
legend("OpenMC results", "Linear fit", "location", "southwest");
set(gca, "FontSize", 12);

print("results/doppler_reactivity_sqrtT.png", "-dpng", "-r300");

% ----------------------------
% Report coefficients
% ----------------------------

disp("Doppler regression summary:");
disp(["  dk/dT (linear fit) = ", num2str(dkdT, "%.4e"), " per K"]);
disp(["  dk/dT (linear fit) = ", num2str(dkdT_pcm, "%.3f"), " pcm/K"]);
disp(["  d(rho)/dT (linear fit) = ", num2str(drhodT_pcm, "%.3f"), " pcm/K"]);
disp(["  d(rho)/d(sqrt(T)) (linear fit) = ", num2str(drho_dsqrtT_pcm, "%.3f"), " pcm/sqrt(K)"]);