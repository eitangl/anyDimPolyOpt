% Produce the numerical results in Example 6.1 of the paper
% Here  u_n = inf_{x \in [-1,1]^n} 5*mean(x)^2 - 4*mean(x.^2)*mean(x) - mean(x)
% and ell_n = inf_{x \in [-1,1]^n} 5*\bar m_n^{(1,1)}(x) - 4\bar m_n^{(2,1)}(x) - \bar m_n^{(1)}(x)
% We apply an SOS relaxation to produce lower bounds for ell_n, and fmincon
% with several initializations to produce upper bounds on u_n and ell_n

clear all, clc
d = 3; % degree of objective
d_sos = 4; % degree of SOS relaxation
n_max = 21; % largest dimension n for u_n and ell_n

% options for fmincon:
opts_fmincon = optimoptions('fmincon','SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','none');
num_inits = 100; % number of initializations for fmincon

% options for mosek (used for SOS):
ops_mosek = sdpsettings('solver','mosek','verbose',0,'debug',0,'cachesolvers',1);

%% produce SOS lower bounds on ell_n 
lower_bounds_sdp = zeros(n_max-d+1,1);

lower = sdpvar(1,1);
for n=d:n_max
    x = sdpvar(n,1);
    
    % form objective polynomial of ell_n
    p = -sum(x)/n;
    for i=1:n-1
        p = p + (5/nchoosek(n,2))*x(i)*sum(x(i+1:end));
    end
    for i=1:n
        p = p - (4/(n*(n-1)))*x(i)^2*sum(x([1:i-1,i+1:n]));
    end
    
    % add SOS multipliers for [-1,1]^n = {1+x>=0, 1-x>=0}
    g = [1+x; 1-x];
    F = [];
    c_vec=[];
    for i=1:length(g)
        [s,c] = polynomial(x,d_sos - 2*ceil(degree(g(i))/2));
        c_vec = [c_vec; c];
        p = p - g(i)*s;
        F = [F; sos(s)];
    end
    
    p = p-lower;
    F = [F;sos(p)];
    solvesos(F, -lower, ops_mosek, [lower;c_vec]);

    lower_bounds_sdp(n-d+1) = value(lower);
end

%% Compute upper bounds on ell_n and u_n using fmincon

lower_bounds = zeros(n_max-d+1,1);
upper_bounds = zeros(n_max-d+1,1);
for n=d:n_max
    % setup box constraints [-1,1]^n
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -ones(n,1);
    ub = ones(n,1);

    min_found_lower = Inf;
    min_found_upper = Inf;
    for t=1:num_inits
        % initialize uniformly on the vertices of the hypercube
        x0 = sign(randn(n,1));

        %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
        [x,fval,exitflag,output] = fmincon(@poly_costgrad_lower, x0,A,b,Aeq,beq,lb,ub,[],opts_fmincon);
        
        if fval < min_found_lower
            min_found_lower = fval;
        end

        [x,fval,exitflag,output] = fmincon(@poly_costgrad_upper, x0,A,b,Aeq,beq,lb,ub,[],opts_fmincon);
        
        if fval < min_found_upper
            min_found_upper = fval;
        end
    end
    lower_bounds(n-d+1) = min_found_lower;
    upper_bounds(n-d+1) = min_found_upper;
end

%% Produce the plots in the paper

lower_bnds_avg = (lower_bounds+lower_bounds_sdp)/2;
lower_bnds_err = (lower_bounds-lower_bounds_sdp)/2;

figure,errorbar(d:n_max, upper_bounds - lower_bnds_avg, lower_bnds_err,'linewidth',2)
hold on
theory_bd = 10*(3*2./(d:n_max)); % theoretical bound on u_n-v_n
plot(d:n_max, theory_bd,'linewidth',2)
%loglog(d:d:n_max, upper_bounds(1:d:end)-lower_bounds(1:d:end),'marker','+','markersize',10,'markeredgecolor','k','linestyle','--','linewidth',1,'color','k')
errorbar(d:d:n_max, upper_bounds(1:d:end)-lower_bnds_avg(1:d:end), lower_bnds_err(1:d:end), 'marker','+','markersize',10,'markeredgecolor','k','linestyle','--','linewidth',1,'color','k')
set(gca, 'XScale','log', 'YScale','log')
legend({'$u_n-\ell_n$','Theor. bd.','$u_{3n}-\ell_{3n}$'},'Interpreter','latex','box','off','location','southwest')
set(gca, 'fontsize',24)
xlabel('$n$','interpreter','latex')
ylabel('$u_n-\ell_n$','interpreter','latex')

figure, plot(d:n_max, upper_bounds, 'x', 'linewidth',2 ,'linestyle','-')
hold on
errorbar(d:n_max, lower_bnds_avg, lower_bnds_err, 'o','linewidth',2,'linestyle', '-')
set(gca, 'fontsize',24)
xlim([d,n_max])
xlabel('$n$','interpreter','latex')
ylabel('Optimal values','interpreter','latex')
legend({'$u_n$', '$\ell_n$'},'interpreter','latex','Location','southeast','box','off')

%% helper funcs

function [f,g] = poly_costgrad_upper(x)
% Calculate objective f
n = length(x);
pwr_means=zeros(2,1);
for i=1:2
    pwr_means(i) = sum(x.^i)/n;
end

f = 5*pwr_means(1)^2 - 4*pwr_means(1)*pwr_means(2) - pwr_means(1);

if nargout > 1 % gradient required
    g = (10*pwr_means(1)/n-4*pwr_means(2)/n-1/n)*ones(n,1) - 8*pwr_means(1)*x/n;
end
end

function [f,g] = poly_costgrad_lower(x)
n = length(x);

f = -sum(x)/n;
for i=1:n-1
    f = f + (5/nchoosek(n,2))*x(i)*sum(x(i+1:end));
end
for i=1:n
    f = f - (4/(n*(n-1)))*x(i)^2*sum(x([1:i-1,i+1:n]));
end

if nargout > 1 % gradient required
    g = (5/nchoosek(n,2))*(sum(x)*ones(n,1) - x) -(4/(n*(n-1)))*(2*sum(x)*x - 3*x.^2 + sum(x.^2)*ones(n,1)) -ones(n,1)/n;
end
end

