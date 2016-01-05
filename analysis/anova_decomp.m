clear all
load('zeta.mat')
countryNames = importdata('country_name.txt');

[d1 d2 d3 d4] = size(zeta);
NN = length(zeta(:));
g1 = repmat([1:d1]', [d2 * d3 * d4, 1]);
g2 = repmat(kron([1:d2]', ones(d1, 1)), [d3 * d4, 1]);
g3 = repmat(kron([1:d3]', ones(d1 * d2, 1)), [d4, 1]);

[p,tbl,stats,terms] = anovan(zeta(:), {g1 g2 g3}, 'varnames',...
    {'importer','exporter','sector'});

coeffs = stats.coeffs;

importer = coeffs(2:26);
exporter = coeffs(27:51);
sector = coeffs(52:74); 

figure()
bar(importer)
title('Importer effect')
set(gca,'XTick', 1:25)
set(gca,'XTickLabel', countryNames)
rotateXLabels(gca(), 60)
figure()
bar(exporter)
title('Exporter effect')
set(gca,'XTick', 1:25)
set(gca,'XTickLabel', countryNames)
rotateXLabels(gca(), 60)
figure()
bar(sector)
title('Sector effect')
set(gca,'XTick', 1:24)

