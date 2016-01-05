function d = correctexpenditureshares(d, parameters)
% There are a number of manipulations to the shares...

numericalZero = parameters.numericalZero;

[nCountries, ~, nSectors, nYears] = size(d);

% There is no data about domestic shares. We calculate them so that
% shares add up to 1.
for n = 1:nCountries
    d(n, n, :, :) = ones(1, 1, nSectors, nYears) - sum(d(n, :, :, :), 2);
end

% (Potentially) We created some negative domestic shares that we have to
% correct for.

% In the first period we clip them from below to 0 and then renormalize the
% shares so that they add up to 1 for a given (destination) country, sector
% time period triplet.
d(d(:, :, :, 1) < 0) = 0;
d(:, :, :, 1) = ...
    d(:, :, :, 1) ./ repmat(sum(d(:, :, :, 1), 2), [1, nCountries, 1]);

% In the later years we approximate the negative domestic share with its
% previous value and then renormalize again.
for t = 2:nYears
    dCurrentYear = d(:, :, :, t);
    dPreviousYear = d(:, :, :, t - 1);
    indexOfNegative = (dCurrentYear < 0);
    dCurrentYear(indexOfNegative) = dPreviousYear(indexOfNegative);
    d(:, :, :, t) = ...
        dCurrentYear ./ repmat(sum(dCurrentYear, 2), [1, nCountries, 1]);
end

% We will take logs and divide by shares, so we set a numerical zero
d(:, :, 1:(nSectors - 1), :) = ...
    max(d(:, :, 1:(nSectors - 1), :), numericalZero);
end