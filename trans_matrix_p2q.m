function [T,partitions] = trans_matrix_p2q(d,k,b1,b2)
%   Produce the transition matrix T converting coefficients of the freely-described cost p_n in basis b1 
%   into coefficients in basis b2 for the freely-symmetrized cost q_k in dimension k (which can be a symbolic variable). 
%   The six bases for symmetric polynomials are
%   e: elementary symmetric function
%   h: complete homogeneous symmetric function
%   m: monomial symmetric function
%   f: forgotten symmetric function
%   s: Schur function
%   p: power-sum symmetric function
%   See also the documentation for: 
%   https://www.mathworks.com/matlabcentral/fileexchange/136299-transition-matrices-between-symmetric-polynomials
%   which is used by the present function.
%
%   Input: 
%   d: degree of the (homogeneous) symmetric polynomials
%   k: dimension for the lower bound (could be symbolic)
%   b1: base 1 for input freely-described polynomial
%   b2: base 2 for output freely-symmetrized polynomial q_k
%   b1 and b2 are either 'e', 'h', 'm', 'f', 's' or 'p'.
%   Output
%   T: transition matrix s.t. if c_p are coefficients of (p_n) in basis b1, 
%       then T*c_p are coefficients of q_k in basis b2
%   partitions: matrix whose rows are partitions on d. The bases b1 and b2
%       are indexed by these partitions, and the order of the coefficients in
%       these bases should be the same as the order of the rows in partitions.
% Eitan Levin, June 6 2025.

T1 = sym(trans_matrixd(d, b1, 'm')); % if v_p are coefficients in b1 basis, then v_p'*T1 are coeffs in m basis
T2 = sym(trans_matrixd(d, b2, 'm'));
R = sym(trans_matrixd(d, 'p','m')); 

partitions = ip_desc(d); % form partitions of d in the right order
lens = sum(partitions > 0, 2); % lengths of these partitions
% compute factorials of partitions
facs = ones(size(lens));
for i=1:size(partitions,1)
    for j=1:lens(i)
        facs(i) = facs(i)*factorial(partitions(i,j));
    end
end

% compute the factoris in Prop. 6.2 of the paper depending on lambda and on
% mu
right_factor = sym(zeros(length(R),1));
left_factor = sym(zeros(length(R),1));
for i=1:length(R)
    right_factor(i) = facs(i)*prod(k-(0:lens(i)-1))/sym(R(i,i));
    left_factor(i) = facs(i)*k^(lens(i));
end

% form the final transition matrix
T = simplify(inv(((left_factor.^(-1)).*R.*right_factor')*T2')*T1', 100);

