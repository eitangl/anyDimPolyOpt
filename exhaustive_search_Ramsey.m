% Produce numerical results for Example 6.7 in the paper. 
% Here u_n = inf_{X simple graph on [n]} t(H;X) + t(H; 1-X) where H=triangle with pendant node,
% and ell_n = ... same but with t_inj...
% We comput ell_n for n=4,...,10 by exhaustive search over simple graphs on
% up to 10 vertices. We obtained a list of these graphs, which the user
% should download and make available in the path, from: https://users.cecs.anu.edu.au/~bdm/data/graphs.html

clear all, clc
A_perms = all_graph_perms([0, [1,0,0]; [1;0;0], ones(3)-eye(3)]); % all permutations of H
k = 4; % number of vertices in H
num_perms = size(A_perms, 3); % size of S_k orbit of H

n_max = 10; % largest n over which to search
lower_bnds = zeros(n_max-k+1,1);
for n = k:n_max
    N = nchoosek(n,2);
    S = nchoosek(1:n, k);
    coeff_H = num_perms*nchoosek(n,k); % normalization for t_inj

    % compute indices to extract from each adjacency matrix to multiply
    % together and form X^{H'} for each H' in S_n orbit of H
    inds = cell(nchoosek(n,k), num_perms); 
    
    x_inds = (1:N)'; % indices of upper triangle
    X_inds = zeros(n);
    X_inds(triu(ones(n),1)==1) = x_inds;
    for j = 1:num_perms
        for s = 1:nchoosek(n,k)
            X = X_inds(S(s,:), S(s,:));
            inds{s,j} = X(triu(ones(k),1)==1 & A_perms(:,:,j)==1);
        end
    end

    % read all graphs on n vertices from file
    graphs = strsplit(fileread(['graph' num2str(n) '.g6']));
    graphs = cellfun(@double, graphs, 'UniformOutput',0);
    if isempty(graphs{end})
        graphs = graphs(1:end-1);
    end
    
    x = zeros(N,1); % variable to hold vectorized entries in upper triangle of adjacency matrices 
    min_found = 1;
    for g = 1:length(graphs)
        % decode the graphs from the file, stored in graph6 format: https://users.cecs.anu.edu.au/~bdm/data/formats.html
        M=ceil(N/6);
        for j = 1:M
            if j < M
                x((j-1)*6+(1:6)) = num2bin(graphs{g}(j+1)-63, 6);
            else
                x_tmp = num2bin(graphs{g}(j+1)-63, 6);
                x((j-1)*6+1:end) = x_tmp(1:N-(j-1)*6);
            end
        end
        
        % form t_inj(H;X) + t_inj(H;1-X)
        m = 0;
        for j = 1:num_perms
            for s = 1:nchoosek(n,k)
                X_curr = x(inds{s,j});
                m = m + prod(X_curr) + prod(1-X_curr);
            end
        end
        m = m/coeff_H;

        if m < min_found
            min_found = m;
        end
        if min_found == 0
            break;
        end
    end
    lower_bnds(n-k+1) = min_found;
end

display('lower bounds:')
lower_bnds'

%% auxiliary functions
function digits = num2bin(n,l)
digits = zeros(l,1);
for j = 1:floor(log2(n))+1
    n = n/2;
    if n ~= floor(n)
        n = floor(n);
        digits(l-j+1) = 1;
    end
end
end