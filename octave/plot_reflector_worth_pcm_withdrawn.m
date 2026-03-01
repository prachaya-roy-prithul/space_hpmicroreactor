% Reflector worth plot

data = dlmread('results/reflector_sweep_withdrawn_clean.csv', ',', 0, 0);

thickness = data(:,2);
k         = data(:,3);
unc       = data(:,4);

rho = (k - 1.0) ./ k;
rho0 = rho(1);

drho = rho - rho0;
worth_pcm = 1.0e5 * drho;

drhodk = 1.0 ./ (k .^ 2);
rho_sd = abs(drhodk) .* unc;
rho0_sd = rho_sd(1);

worth_pcm_sd = 1.0e5 * sqrt(rho_sd.^2 + rho0_sd.^2);

figure('color','w');

h = errorbar(thickness, worth_pcm, worth_pcm_sd, 'o-');
set(h, 'linewidth', 1.5);
set(h, 'markersize', 7);

xlim([0 12]);
ylim([0 42000]);

xlabel('Reflector thickness (cm)', 'fontsize', 12);
ylabel('Reflector worth, \Delta\rho (pcm)', 'fontsize', 12);

title({'Reflector worth vs thickness', ...
       '(T_{fuel} = 600 K, control rod withdrawn)'}, ...
       'fontsize', 13);

grid on;
set(gca, 'gridalpha', 0.15);

box on;

print('-dpng', '-r300', 'results/reflector_worth_pcm_withdrawn.png');

fprintf('Saved reflector_worth_pcm_withdrawn figure to results/reflector_worth_pcm_withdrawn.png\n');