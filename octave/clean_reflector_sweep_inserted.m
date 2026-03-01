% Clean reflector sweep CSV (inserted) into numeric table for dlmread workflows.
% Output columns:
%   [Tfuel_K, reflector_thickness_cm, keff, unc, leak_frac, leak_frac_unc]

in_path  = 'results/reflector_sweep_inserted.csv';
out_path = 'results/reflector_sweep_inserted_clean.csv';

fid = fopen(in_path, 'r');
if fid < 0
  error(['Unable to open file: ', in_path]);
end

header = fgetl(fid); %#ok<NASGU>

C = textscan(fid, '%s %f %f %f %f %f %f', 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);

if isempty(C) || isempty(C{2})
  error('No numeric data parsed from inserted CSV.');
end

num = C{2};

Tfuel = num(:,1);
refl  = num(:,2);
keff  = num(:,3);
unc   = num(:,4);
leak  = num(:,5);
leaku = num(:,6);

M = [Tfuel, refl, keff, unc, leak, leaku];

valid = isfinite(M(:,1)) & isfinite(M(:,2)) & isfinite(M(:,3)) & isfinite(M(:,4)) & isfinite(M(:,5)) & isfinite(M(:,6));
M = M(valid, :);

if rows(M) == 0
  error('All rows removed during cleaning (no finite rows).');
end

% Sort by thickness ascending, and keep the LAST occurrence of each thickness
[~, order] = sort(M(:,2));
M = M(order, :);

kept = true(rows(M), 1);
for r = 1:(rows(M)-1)
  if abs(M(r,2) - M(r+1,2)) < 1.0e-12
    kept(r) = false;
  end
end
M = M(kept, :);

dlmwrite(out_path, M, ',', 'precision', 12);

fprintf('Wrote %s (%d rows)\n', out_path, rows(M));