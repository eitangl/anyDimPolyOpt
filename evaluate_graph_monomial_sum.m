function m = evaluate_graph_monomial_sum(A_perms, X)
% Given all permutations of the adjacency matrix and an input adjacency
% matrix on at least as many vertices, compute the graph monomial sums
% m^[H](X).
% Inputs:
% A_perms: tensor of size k x k x |Aut(H)| of all permutations of the
%   adjacency matrix of a graph H on k vertices, obtained using
%   all_graph_perms, for instance.
% X: symmetric n x n matrix with n>=k
% Outputs:
% m: value of m_n^[H](X)
% Eitan Levin, June 6 2024.

k = size(A_perms, 1);
n = size(X,1);
S = nchoosek(1:n, k);
m = 0;
vec = @(x) x(:);
for i = 1:size(S,1)
    for j = 1:size(A_perms,3)
        X_curr = triu(X(S(i,:),S(i,:))) + tril(ones(k),-1);
        A_curr = A_perms(:,:,j);
        m = m + prod(vec(X_curr.^(A_curr)));
    end
end