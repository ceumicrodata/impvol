function e = get_eq_trade(model)

load(['models/', model, '/equilibrium'])
load(['models/', model, '/alg_inputs'])
load(['models/', model, '/data_rgdp_and_volatility.mat'], 'd', 'va')

e.data_d = d;
clear d

% load import shares
d = equilibrium.d;
d_imp = d;
for t = 1:36
    for j = 1:24
        a = d_imp(:, :, j, t);
        a(logical(eye(25))) = 0;
        d_imp(:, :, j, t) = a;
    end
end

% load equilibrium quantities
L_njt = equilibrium.L_njt;
w_njt = equilibrium.w_njt;
% P_njt = equilibrium.P_njt;
% L_nt = equilibrium.L_nt;
% w_nt = equilibrium.w_nt;
% P_nt = equilibrium.P_nt;

% load i-o links
gamma = baseline.gammas;

% load alphas and betas
alpha = baseline.alpha;
beta = baseline.beta;

% load trade balance
S = baseline.trade_balance;

% compute expenditure
PC_nt = squeeze(sum(L_njt .* w_njt, 2)) - S;

R_njt = L_njt .* w_njt ./ permute(repmat(beta, [1, 25, 36]), [2, 1, 3]);

E_njt = permute(repmat(alpha, [1, 1, 25]), [3, 1, 2]) .* ...
        permute(repmat(PC_nt, [1, 1, 24]), [1, 3, 2]);
    
for n = 1:25
    for j = 1:24
        for t = 1:36
            E_njt(n, j, t) = E_njt(n, j, t) + gamma(j, :, 1) * R_njt(n, :, t)';
        end
    end
end


va_njt = L_njt .* w_njt;
va_nt = squeeze(sum(va_njt, 2));
e.va_t = sum(va_nt, 1)';

e.I_njt = squeeze(sum(d_imp, 2)) .* E_njt;
e.I_nt = squeeze(sum(e.I_njt, 2));
e.I_t = sum(e.I_nt, 1)';

% trade_share = im_t ./ va_t;

e.alpha = alpha;
e.beta = beta;
e.gamma = gamma;
e.S = S;

e.va_njt = va_njt;
e.va_nt = va_nt;
e.E_njt = E_njt;
e.E_nt = squeeze(sum(E_njt, 2));

e.d = d;
e.d_imp = d_imp;

e.omd_njt = squeeze(sum(e.d_imp, 2));
e.omd_nt = squeeze(sum(e.omd_njt, 2));
e.omd_t = sum(e.omd_nt, 1)';
e.names = baseline.country_names;

e.kappa = baseline.kappa;


dd_imp = e.data_d;
for t = 1:36
    for j = 1:24
        a = dd_imp(:, :, j, t);
        a(logical(eye(25))) = 0;
        dd_imp(:, :, j, t) = a;
    end
end

e.dd_imp = dd_imp;

e.omdd_njt = squeeze(sum(e.dd_imp, 2));
e.omdd_nt = squeeze(sum(e.omdd_njt, 2));
e.omdd_t = sum(e.omdd_nt, 1)';

%% Data trade flows
% compute expenditure
PC_nt = squeeze(sum(va, 2)) - S;

R_njt = va ./ permute(repmat(beta, [1, 25, 36]), [2, 1, 3]);

E_njt = permute(repmat(alpha, [1, 1, 25]), [3, 1, 2]) .* ...
        permute(repmat(PC_nt, [1, 1, 24]), [1, 3, 2]);
    
for n = 1:25
    for j = 1:24
        for t = 1:36
            E_njt(n, j, t) = E_njt(n, j, t) + gamma(j, :, 1) * R_njt(n, :, t)';
        end
    end
end

e.dva_nt = squeeze(sum(va, 2));

e.dI_njt = squeeze(sum(dd_imp, 2)) .* E_njt;
e.dI_nt = squeeze(sum(e.dI_njt, 2));
e.dI_t = sum(e.dI_nt, 1)';
