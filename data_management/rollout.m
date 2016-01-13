function out = rollout(change, level)

[n_countries, n_sectors, n_years] = size(level);

out = zeros(n_countries, n_sectors, n_years);
out(:, :, 1) = level(:, :, 1);

for t = 2:n_years
    out(:, :, t) = out(:, :, t - 1) .* exp(change(:, :, t - 1));
end

end