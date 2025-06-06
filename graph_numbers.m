% Produce numerical results for Example 6.10 in the paper. 
% Here u_n = inf_{X \in S^n, X(:) >= 0, sum(X(:))=1} inj(P_3; X)-inj(K_2+K_2, X),
% and ell_n = inf_{same constraints} (n^3/(n(n-1)(n-2)))*[ inj(P_3;X) - (n+1)/(n-3)inj(K_2+K_2; X) ]
% We produce SOS lower bounds on ell_n as well as upper bounds on ell_n and u_n using fmincon
% The code below can perform these computations for any linear combination
% of graph monomial sums m^[H] = inj(H; .) / |Aut(H)|
clear all, clc

A_in = {[0,1,0;1,0,1;0,1,0], [0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0]}; % adjacency matrices of the two graphs above
c_in = @(k) [2*k^3/(k*(k-1)*(k-2)), -8*k^3*(k+1)/(k*(k-1)*(k-2)*(k-3))]; % coefficients of corresponding monomial sums in ell_n
c_in_upper = [2, -8]; % coefficients of corresponding monomial sums in u_n

A_in = cellfun(@all_graph_perms, A_in, 'UniformOutput',false); % form orbits of graphs above
num_verts = max(cellfun(@(x)size(x,1), A_in)); % get their number of vertices
num_perms = max(cellfun(@(x)size(x,3), A_in)); % get the size of their orbits

k=max(num_verts); % max number of vertices
n_max = 9; % max dimension n to consider
d_sos = 4; % degree of SOS relaxation
% options for mosek (used for SOS)
ops = sdpsettings('solver','mosek','verbose',2,'debug',0,'cachesolvers',1);
% options for fmincon:
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','none');
num_inits = 5; % number of fmincon initializataions

%% produce lower bounds using SOS
lower_bounds_sdp = zeros(n_max-k+1,1);

lower = sdpvar(1,1);
for n=k:n_max
    c_curr = c_in(n);
    X = sdpvar(n);
    X(eye(n)==1)=0; % we can restrict to zero-diagonal matrices since graphs in A_in are simple
    x_triu = X(triu(ones(n),1)==1);
    
    % form the linear combination of graph monomial sums
    p = 0;
    for i=1:length(A_in)
        m_H = evaluate_graph_monomial_sum(A_in{i}, X);
        p = p + c_curr(i)*m_H;
    end
    
    % add SOS multipliers for X>=0 constraint
    g = [x_triu];
    F = [];
    c_vec = [];
    mons = monolist(x_triu, d_sos-2);
    for i=1:length(g)
        c = sdpvar(length(mons),1);
        s = c'*mons;
        p = p - g(i)*s;
        F = [F; sos(s)];
        c_vec = [c_vec; c];
    end
    p = replace(p, x_triu(1), 1-sum(x_triu(2:end))); % enforce sum(X(:))=1 constraint
    p = p - lower;
    F = [F; sos(p)];

    solvesos(F, -lower, ops, [lower;c_vec]);

    lower_bounds_sdp(n-k+1) = value(lower);
end

%% produce upper bounds using fmincon

lower_bounds = zeros(n_max-k+1,1);
upper_bounds = zeros(n_max-k+1,1);
for n = k:n_max
    % setup simplex constraint
    N = nchoosek(n,2);
    A = [];
    b = [];
    Aeq = ones(1,N);
    beq = 1;
    lb = zeros(N,1);
    ub = ones(N,1);

    min_found_lower = Inf;
    min_found_upper = Inf;
    for t=1:num_inits
        % initialize uniformly on simplex
        x0 = exprnd(1,[N,1]);
        x0 = x0/sum(x0);
        %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
        [x,fval,exitflag,output] = fmincon(@(x) poly_cost(x, A_in, c_in(n), n), x0,A,b,Aeq,beq,lb,ub,[],options);
        
        if fval < min_found_lower
            min_found_lower = fval;
        end

        [x,fval,exitflag,output] = fmincon(@(x) poly_cost(x, A_in, c_in_upper, n), x0,A,b,Aeq,beq,lb,ub,[],options);
        
        if fval < min_found_upper
            min_found_upper = fval;
        end
    end
    lower_bounds(n-k+1) = min_found_lower;
    upper_bounds(n-k+1) = min_found_upper;
end

%% Display results

display('Table with rows n, upper bound on u_n, upper bound on ell_n, lower bound on ell_n:')
[(d:n_max)',upper_bounds,lower_bounds,lower_bounds_sdp]'


%% Helper function
function [p,g] = poly_cost(x, A_in, c_in, n)
% Calculate objective f
X = zeros(n);
X(triu(ones(n),1)==1) = x;
X = X+X'; 
% X(eye(n)==1) = X(eye(n)==1)/2;

p = 0;
for i=1:length(A_in)
    m_H = evaluate_graph_monomial_sum(A_in{i}, X);
    p = p + c_in(i)*m_H;
end
if nargout > 1 % gradient required
    sx = sum(X,2);
    g = (c_in(1) - c_in(2))*(sx*ones(1,n) + ones(n,1)*sx')...
        + c_in(2)*(sum(x)) - (2*c_in(1)-c_in(2))*X;
    g = g(triu(ones(n),1)==1);
end
end