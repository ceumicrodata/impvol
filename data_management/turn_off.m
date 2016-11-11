function d = turn_off(n, d)

[~, ~, n_sectors, n_years] = size(d);

for t = 1:n_years
    for j = 1:n_sectors
        d(n, :, j, t) = 0;
        d(:, n, j, t) = 0;
%         d(n, n, j, t) = 1;
    end
end