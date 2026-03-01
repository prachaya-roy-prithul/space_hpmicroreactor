% Reflector sweep CSV cleaning utility

in_csv  = 'results/reflector_sweep_withdrawn.csv';
out_csv = 'results/reflector_sweep_withdrawn_clean.csv';

fid = fopen(in_csv);

% Read header
fgetl(fid);

% Read data: skip string column, read numeric columns
C = textscan(fid, '%s %f %f %f %f %f %f', 'Delimiter', ',');

fclose(fid);

Tfuel     = C{2};
thickness = C{3};
keff      = C{4};
unc       = C{5};
leak      = C{6};
leak_unc  = C{7};

% Remove NaN leakage rows
mask = ~isnan(leak);
Tfuel     = Tfuel(mask);
thickness = thickness(mask);
keff      = keff(mask);
unc       = unc(mask);
leak      = leak(mask);
leak_unc  = leak_unc(mask);

% Remove duplicate thickness values (keep first)
[thickness_unique, idx] = unique(thickness, 'stable');

Tfuel     = Tfuel(idx);
keff      = keff(idx);
unc       = unc(idx);
leak      = leak(idx);
leak_unc  = leak_unc(idx);

clean_data = [Tfuel thickness_unique keff unc leak leak_unc];

dlmwrite(out_csv, clean_data, 'delimiter', ',', 'precision', 10);

fprintf('Cleaned sweep written to %s\n', out_csv);