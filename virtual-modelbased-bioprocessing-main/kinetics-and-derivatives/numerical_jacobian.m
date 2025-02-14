%% Function to evaluate a numerical Jacobian by finite differences
%   author        : Mirko Pasquini 
%
%  The function takes as input a vector function handler, a point in which
%  the Jacobian should be evaluated and an epsilon for the numerical
%  differentiation. It will return a Jacobian evaluated by numerical
%  differences, i.e. J(i,j) = (f_i(u0+e_j*h_j)-f(u0-e_j*h_j))/(2*h_j) where
%  e_j is the j-th column of the identity matrix and h_j is the j-th
%  element of the vector h (see definition below)
%
%  J_num = NUMERICAL_JACOBIAN(fun,u0,h)
%
%  @inputs:
%       fun : handle of a vector function 
%       u0 : point in which the Jacobian should be evaluated
%       h : epsilons vector for numerical differentiation. Should be
%           the same size of the vector u0
%
%  @outputs:
%       J_num : Jacobian evaluated by finite differences

function J_num = numerical_jacobian(fun,u0,h)
    f0 = fun(u0);
    nf = length(f0);
    nu = length(u0);
    J_num = zeros(nf,nu);
    for j = 1 : nu
        hj = h(j);
        uj_plus = u0;
        uj_minus = u0;
        uj_plus(j) = u0(j)+hj;
        uj_minus(j) = u0(j)-hj;
        J_num(:,j) = (fun(uj_plus)-fun(uj_minus))/(2*hj);
    end
end