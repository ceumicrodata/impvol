function out = rollout(change, level)

[nCountries, nSectors, nYears] = size(level);

out = zeros(nCountries, nSectors, nYears);
out(:, :, 1) = level(:, :, 1);

for t = 2:nYears
    out(:, :, t) = out(:, :, t - 1) .* exp(change(:, :, t - 1));
end

end