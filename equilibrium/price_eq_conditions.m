function [Pmtheta_nj_new] = price_eq_conditions(Pmtheta_nj, D, t)

global B xi gammas theta

[N, ~, J] = size(D);

aa = reshape(repmat(prod(kron(Pmtheta_nj', ones(1, J)) .^ repmat(gammas(:, :, t), [1, N])), [N, 1]), [J * N, N]);
bb = reshape(sum(aa .* reshape(permute(D, [1, 3, 2]), [N * J, N]), 2), [N, J]);
Pmtheta_nj_new = xi^(-theta) * bsxfun(@times, bb, B(:, t)'.^(-theta));


% Pmtheta_nj_new2 = zeros(N, J);
% 
% for n = 1:N
%     for j = 1:J
%         
%         sum_over_n = 0;
%         for i = 1:N
%             
%             prod_P = 1;
%             
%             for k = 1:J
%                 prod_P = prod_P * Pmtheta_nj(i, k)^gammas(k, j);
%             end % for k
%             
%             sum_over_n = sum_over_n + D(n, i, j) * prod_P;
%         end % for i
%         
%         Pmtheta_nj_new2(n, j) = (xi * B(j))^(- theta) * sum_over_n;% - Pmtheta_nj(n, j);
%             
%     end % for j
% end % for n
% 
% maxdif = max(abs(Pmtheta_nj_new(:) - Pmtheta_nj_new2(:)));


% jacobian

% for n = 1:N
%     for j = 1:J
%         
%         row = (j - 1) * N + n;
%         
%         for m = 1:N
%             for s = 1:J
%                 col = (s - 1) * N + m;
%                 
%                 prod_P = 1;
%                 for k = 1:J
%                     prod_P = prod_P * Pmtheta_nj(m, k)^gammas(k, j);
%                 end % for k
%                 
%                 jacobian(row, col) = (xi * B(j))^(- theta) * D(n, m, j) * gammas(s, j) * prod_P / Pmtheta_nj(m, s);
%                 
%             end % for s
%         end % for i
%     end % for j
% end % for n
% 
% jacobian = jacobian - eye(N * J);