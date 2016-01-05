function out = price_contraction(P, D, t)
% also uses globals: alpha, beta, K
% N nonlinear equations of N unknowns in the form of F(P)

% Inputs:
% P: N x 1 vector of prices
% alpha, beta, : J x 1 vectors of alpha and beta values
% K: constant
% D: N x N x J array of coefficients

global alpha beta K

N = length(P);
[J , ~] = size(alpha);

% Output:
% out is an N x 1 vector

% preallocate memory for the N x J sums: (sum(D_ij * P_i^(beta_j - 1)))^(- alpha_j)
sumDPba = zeros(N, J);

for j = 1:J
    sumDPba(:, j) = (D(:, :, j) * P.^(beta(j) - 1)).^(- alpha(j, t));
end

out = K(t) * prod(sumDPba, 2);
end