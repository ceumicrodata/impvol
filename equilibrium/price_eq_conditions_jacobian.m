function fn_value = price_eq_conditions_jacobian(Pmtheta_nj, D)

global B xi gammas theta

[N, ~, J] = size(D);

% benchmark. very inefficient. time: cca. 0.04 sec
fn_value = zeros(N, J);

for n = 1:N
    for j = 1:J
        
        sum_over_n = 0;
        for i = 1:N
            
            prod_P = 1;
            
            for k = 1:J
                prod_P = prod_P * Pmtheta_nj(i, k)^gammas(k, j);
            end % for k
            
            sum_over_n = sum_over_n + D(n, i, j) * prod_P;
        end % for i
        
        fn_value(n, j) = (xi * B(j))^(- theta) * sum_over_n - Pmtheta_nj(n, j);
            
    end % for j
end % for n

% benchmark. very inefficient. time: cca. 1.25 sec
% jacobian = zeros(N * J);

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