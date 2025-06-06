% Produce numerical results for Example 6.4 in the paper.
% Here u_n = \inf_{|x|_1 <= 1} p_n(x) = 4s_n^{(4)}(x) - (139/20)s_n^{(3,1)}(x) + 4s_n^{(2,2)}(x) - 5s_n^{(2,1,1)}(x) + 4s_n^{(1,1,1,1)}(x)
% and ell_n = \inf_{|x|_1 <= 1} q_n(x) where the coefficients of q_n in the power sum basis are generated computationally below.

d = 4; % degree of cost
n_sym = sym("n","integer"); % symbolic variable for coefficients of lower bound

% generate matrix converting from coefficients of p_n(x) in power sum basis
% to coefficients of q_{n_sym}(x) in power sum basis, where both
% coefficients are listed in the order of the rows of partitions
[T,partitions] = trans_matrix_p2q(d,n_sym,'p','p'); 
c_in = [4,-sym(139/20),4,-5,4]'; % coefficients of p_n
c_out = T*c_in; % coefficients of q_n, as a function of n

d_sos = 6; % degree of SOS relaxation
n_max = 8; % max dimension for lower bound
% options for mosek (used for SOS):
ops = sdpsettings('solver','mosek','verbose',1,'debug',0,'cachesolvers',1);
% options for fmincon:
opts_fmincon = optimoptions('fmincon','SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','none');
num_inits = 10; % number of initializations for fmincon

%% produce lower bounds using SOS
lower_bounds_sdp = zeros(n_max-d+1,1);

lower = sdpvar(1,1);
for n=d:n_max
    c_curr = double(subs(c_out, n_sym, n)); % q_n coefficients for current n

    x = sdpvar(n,1);
    y = sdpvar(n,1); % extra variable for ell_1 ball constraint

    % form power sums
    pwr_sums=[];
    for i=1:d
        pwr_sums = [pwr_sums;sum(x.^i)];
    end
    
    % form q_n
    p = 0;
    for i=1:length(partitions)
        curr_term=1;
        len_curr = sum(partitions(i,:)>0);
        for j=1:len_curr
            curr_term = curr_term*pwr_sums(partitions(i,j));
        end
        p = p + c_curr(i)*curr_term;
    end

    % add SOS multipliers for y-x>=0 and y+x>=0 constraints
    g = [y-x;y+x];
    F = [];
    mons = monolist([x;y], d_sos-2);
    c_vec=[];
    for i=1:length(g)
        c = sdpvar(length(mons),1);
        s = c'*mons;
        c_vec = [c_vec; c];

        p = p - replace(g(i)*s, y(n), 1-sum(y(1:n-1))); % enforce sum(y)=1 constraint
        F = [F; sos(s)];
    end

    p = p-lower;
    F = [F;sos(p)];
    solvesos(F, -lower, ops, [lower;c_vec]);

    lower_bounds_sdp(n-d+1) = value(lower);
end

%% Produce upper bounds on ell_n using fmincon
lower_bounds_fmincon = zeros(n_max-d+1,1);
for n=d:n_max
    A = [eye(n), -eye(n); -eye(n), -eye(n);zeros(1,n), ones(1,n)];
    b = [zeros(2*n,1);1];
    Aeq = [];
    beq = [];
    lb = -ones(2*n,1);
    ub = ones(2*n,1);

    min_found = Inf;
    for t=1:num_inits
        % initialize randomly on the face of the ell_1 ball that seems to
        % achieve the minimum for this objective
        x0 = exprnd(1,[n,1]);
        x0 = x0/sum(x0);
        x0 = x0.*[ones(floor(n/2),1);-ones(n-floor(n/2),1)];
        x0 = [x0; abs(x0)];

        %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
        [x,fval,exitflag,output] = fmincon(@(x) poly_costgrad_ext(x,d,partitions,double(subs(c_out, n_sym, n))),...
            x0,A,b,Aeq,beq,lb,ub,[],opts_fmincon);
        
        if fval < min_found
            min_found = fval;
        end
    end
    lower_bounds_fmincon(n-d+1) = min_found;
end

%% Plot
display('Table with rows n, SOS lower bound, fmincon upper bound:')
[(d:n_max)',lower_bounds_sdp,lower_bounds_fmincon]'

%% helper function
function [f,g] = poly_costgrad_ext(x,d,partitions,coeffs)
n = length(x)/2;
x = x(1:n);
% Calculate objective f
pwr_sums=zeros(d,1);
for i=1:d
    pwr_sums(i) = sum(x.^i);
end

p = 0;
curr_term=ones(length(partitions),1);
for i=1:length(partitions)
    len_curr = sum(partitions(i,:)>0);
    for j=1:len_curr
        curr_term(i) = curr_term(i)*pwr_sums(partitions(i,j));
    end
    p = p + coeffs(i)*curr_term(i);
end
f = p;

if nargout > 1 % gradient required
    g = zeros(size(x));
    for i=1:length(partitions)
        len_curr = sum(partitions(i,:)>0);
        for j=1:len_curr
            g = g + coeffs(i)*(partitions(i,j)*curr_term(i)/pwr_sums(partitions(i,j)))*x.^(partitions(i,j)-1);
        end
    end
    g = [g;zeros(n,1)];
end
end

