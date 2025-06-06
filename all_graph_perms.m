function A_perms = all_graph_perms(A)
% form all permutations of an adjacency matrix A
P = perms(1:length(A));

A_perms = zeros(size(P,1), numel(A));
vec = @(x)x(:);
for i=1:size(P,1)
    A_perms(i,:) = vec(A(P(i,:),P(i,:)))';
end

A_perms = unique(A_perms, 'rows', 'sorted');
A_perms = reshape(A_perms', [size(A,1), size(A,2),size(A_perms,1)]);