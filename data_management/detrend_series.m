function [S_trend, S_cycle] = detrend_series(S, weights)

if ndims(S) == 3
    [N, J, T] = size(S);

    % reshape series to (N * J) x (T)
    S = reshape(S, [N * J, T]);
    
elseif ismatrix(S) == 1
    [N, T] = size(S);
    J = 1;
    
else
    fprintf('Incorrect input in detrend_series.m \n')
    error('Series cannot be detrended.')
end %if

% S_trend = repmat(mean(S, 2), [1, T]);
S_trend = zeros([N * J, T]);
S_cycle = zeros([N * J, T]);

% get maximal lag used in the moving average
K = length(weights) - 1;
weightsVector = [flipud(weights); weights(2:end)];

% augment series by putting in K zeros at edges
%S_front = S(:, 2:(K + 1));
%S_end = S(:, 2:(K + 1))
S_augmented = zeros(N * J, T + 2 * K);
S_augmented(:, (K + 1):(K + T)) = S;

% fill with edges
% S_head_filler = S(:, 1);
% S_tail_filler = S(:, T);

% fill with mean of (K + 1)-window
% S_head_filler = mean(S(:, 1:(K)), 2);
% S_tail_filler = mean(S(:, (T - K + 1):T), 2);

% fit edges
% S_head_filler = S(:, 1:(K + 1)) * (weights / sum(weights));
% S_tail_filler = S(:, (T - K):T) * flipud(weights / sum(weights));

% S_augmented(:, 1:K) = repmat(S_head_filler, [1, K]);
% S_augmented(:, (T + K + 1):(T + 2 * K)) = repmat(S_tail_filler, [1, K]);

% reflect head/tail
S_head_filler = ...
    repmat(S(:, 1), [1, K]) + ...
    fliplr(repmat(S(:, 1), [1, K]) - S(:, 2:(K + 1)));
S_tail_filler = ...
    repmat(S(:, T), [1, K]) + ... 
    fliplr(repmat(S(:, T), [1, K]) - S(:, (T - K):(T - 1)));

S_augmented(:, 1:K) = S_head_filler;
S_augmented(:, (T + K + 1):(T + 2 * K)) = S_tail_filler;



for t = 1:T
    S_cycle(:, t) = S_augmented(:, t:(t + 2 * K)) * weightsVector;
    S_trend(:, t) = S_augmented(:, t + K) - S_cycle(:, t);
end

if J > 1
    S_trend = reshape(S_trend, [N, J, T]);
    S_cycle = reshape(S_cycle, [N, J, T]);
end %if


end