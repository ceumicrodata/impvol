
clear all

%% Baseline
m0 = 'table_th_4_lac_inf_baseline';
e0 = get_eq_trade(m0);

%% Kappa 1972 counterfactual 
m1 = 'table_th_4_lac_inf_kappa1972';
e1 = get_eq_trade(m1);

%% check 1 - S/Y
a = e0.S ./ e0.va_nt;
figure(200000)
hist(1 - a(:), 50)
title('Distribution of (1 - S/Y)')

%% check (1-beta)/beta
gamma = e0.gamma(:, :, 1);
beta = e0.beta;

figure(1)
bar((1-beta)./beta)
title('(1-beta)/beta')
xlabel('Sector')


%% Unweighted d
figure(2)
plot(1972:2007, [e0.omd_t, e1.omd_t])
title('Unweighted d')
grid on
xlabel('year')
legend('Baseline', '1972 costs', 'Location', 'NW')



%% Weighted d: a
alpha = e0.alpha;

e0.omd_a = squeeze(sum(bsxfun(@times, e0.omd_njt, permute(alpha, [3, 1, 2])), 2));
e1.omd_a = squeeze(sum(bsxfun(@times, e1.omd_njt, permute(alpha, [3, 1, 2])), 2));

figure(3)
plot(1972:2007, [mean(e0.omd_a); mean(e1.omd_a)]')
title('d weighted by alpha')
grid on
xlabel('year')
legend('Baseline', '1972 costs', 'Location', 'NW')



%% Weighted d: a + 2
% alpha = e0.alpha;

e0.omd_a2 = squeeze(sum(bsxfun(@times, e0.omd_njt, permute(alpha + 2, [3, 1, 2])), 2));
e1.omd_a2 = squeeze(sum(bsxfun(@times, e1.omd_njt, permute(alpha + 2, [3, 1, 2])), 2));

figure(4)
plot(1972:2007, [mean(e0.omd_a2); mean(e1.omd_a2)]')
title('d weighted by (alpha + 2)')
grid on
xlabel('year')
legend('Baseline', '1972 costs', 'Location', 'NW')



%% Weighted d: a + (1 - b)/b
w = bsxfun(@plus, alpha, (1 - beta) ./ beta);

e0.omd_ab = squeeze(sum(bsxfun(@times, e0.omd_njt, permute(w, [3, 1, 2])), 2));
e1.omd_ab = squeeze(sum(bsxfun(@times, e1.omd_njt, permute(w, [3, 1, 2])), 2));

figure(5)
plot(1972:2007, [sum(e0.omd_ab); sum(e1.omd_ab)]')
title('d weighted by (alpha + (1-b)/b)')
grid on
xlabel('year')
legend('Baseline', '1972 costs', 'Location', 'NW')



%% Actual trade volumes
e0.share_nt = e0.I_nt ./ e0.va_nt;
e1.share_nt = e1.I_nt ./ e1.va_nt;

figure()
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.share_nt(n, :); e1.share_nt(n, :)]')
    title(e0.names{n})
end


tr0 = squeeze(sum(sum(e0.omd_njt .* e0.va_njt, 1), 2)) ./ e0.va_t;
tr1 = squeeze(sum(sum(e1.omd_njt .* e1.va_njt, 1), 2)) ./ e1.va_t;

plot([tr0, tr1])

figure(1)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.va_nt(n, :); e1.va_nt(n, :)]')
    title(e0.names{n})
end

figure(2)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.E_nt(n, :); e1.E_nt(n, :)]')
    title(e0.names{n})
end


figure(3)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e1.E_nt(n, :) ./ e1.va_nt(n, :)]')
    title(e0.names{n})
end


%%

figure(1)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.omd_nt(n, :); e1.omd_nt(n, :)]')
    title(e0.names{n})
end

figure()
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.omd_a(n, :); e1.omd_a(n, :)]')
    title(e0.names{n})
end


figure(2)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.omd_a(n, :); e1.omd_a(n, :)]')
    title(e0.names{n})
end

figure(3)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.omd_ab(n, :); e1.omd_ab(n, :)]')
    title(e0.names{n})
end



figure(4)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.I_nt(n, :); e1.I_nt(n, :)]')
    title(e0.names{n})
end




figure(5)
plot(1972:2007, [e0.I_t, e1.I_t])
title('Total trade volume')
xlabel('year')



%% Kappas

ch = squeeze(e0.kappa(5, :, :, :));


k5 = kappa(:, :, 1:23, :);
aa = squeeze(mean(mean(k5, 2), 3));
plot(1972:2007, )
title('Mean trade cost, Belgium and Luxembourg')

figure(1)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, aa(n,:))
    title(e0.names{n})
end

%%
n = 1;
m = 8;
j = 12;
% k1 = squeeze(kappa(n, m, j, :));
% k2 = squeeze(e0.kappa(n, m, j, :));
k2 = squeeze(e0.d(n, m, j, :));
figure(1)
% plot([k1, k2])
plot(k2)
%%

k5 = kappa(5, :, :, :);
plot(squeeze(mean(mean(k5, 2), 3)))


k36 = e1.kappa(:, 5, :, 36);
mean(k36(:))

plot(squeeze(e0.kappa(1, 5, 1, :)))

n_countries = 25;
n_sectors = 24;
n_years = 36;
i_services = 24;
numerical_zero = 1e-12;

kappa = zeros(n_countries, n_countries, n_sectors, n_years);
theta = 4;
d = e1.d;
for j = 1:(n_sectors - 1)
    for t = 1:n_years
        kappa(:, :, j, t) = ...
            ((d(:, :, j, t) .* d(:, :, j, t)') ./ ...
             (diag(d(:, :, j, t)) * diag(d(:, :, j, t))')).^(1 / (2 * theta));
         
         % normalize
%          kappa(:, :, j, t) = ...
%              kappa(:, :, j, t) ./ ...
%              repmat(max(kappa(:, :, j, t), [], 2), [1, n_countries]);
%          
%          for n = 1:n_countries
%              kappa(n, n, j, t) = 1;
%          end % n
    end % t
end % j


kappa(:, :, i_services, :) = repmat(eye(n_countries), [1, 1, 1, n_years]);

% clip values to the [numerical_zero, 1] interval
kappa(:, :, 1:(n_sectors - 1), :) = ...
    max(kappa(:, :, 1:(n_sectors - 1), :), numerical_zero);
kappa = min(kappa, 1);


aa = e1.kappa - kappa;

mean(aa(:))
std(aa(:))

squeeze(e1.kappa(15, 1:5, 1:5, 1))
squeeze(kappa(15, 1:5, 1:5, 1))


%% Value added (data)

figure(1)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, e0.va_nt(n, :))
    title(e0.names{n})
end


%% Data d

e0.dva_t = sum(e0.dva_nt, 1)';
figure()
plot(1972:2007,[e0.dI_t ./ e0.dva_t, e0.I_t ./ e0.va_t, e1.I_t ./ e1.va_t])
title('Trade share')
grid on
xlabel('year')
legend('Data', 'Baseline', '1972', 'Location', 'NW')

e0.dI_nt(12, :) = [];
e0.dva_nt(12, :) = [];
e0.dI_t = sum(e0.dI_nt, 1)';
e0.dva_t = sum(e0.dva_nt, 1)';

e0.I_nt(12, :) = [];
e0.va_nt(12, :) = [];
e0.I_t = sum(e0.I_nt, 1)';
e0.va_t = sum(e0.va_nt, 1)';

e1.I_nt(12, :) = [];
e1.va_nt(12, :) = [];
e1.I_t = sum(e1.I_nt, 1)';
e1.va_t = sum(e1.va_nt, 1)';


%%
figure(1)
for n = 1:25
    subplot(5, 5, n)
    plot(1972:2007, [e0.dI_nt(n, :) ./ e0.dva_nt(n, :);...
                     e0.I_nt(n, :) ./ e0.va_nt(n, :);...
                     e1.I_nt(n, :) ./ e1.va_nt(n, :)]')
    title(e0.names{n})
    if n == 1
        legend({'Data', 'Baseline', '1972 costs'}, 'Position',[0.07,0.885,0,0])
%         legend('Position',[0.5,0.5,0.5,0])
%         legend('boxoff')
    end
end




