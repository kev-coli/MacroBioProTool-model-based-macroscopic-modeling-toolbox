%  function [x_opt,cost_opt] = quadprog_mosek_as_second_order_programming(K,b,Aineq,bineq,Aeq,beq,xlb,xub)
%
%  This code solves the quadratic programming 
%         min ||K*x-b||^2
%     subject to Aeq*x=beq
%                Aineq*x <= bineq
%                xlb <= x <= xub
%  as a second order programming with MOSEK of the following form
%         min t
%     subject to Aeq*x=beq
%                Aineq*x <= bineq
%                xlb <= x <= xub
%                t >= ||K*x-b||

function [x_opt,cost_opt] = quadprog_mosek_as_second_order_programming(K,b,Aineq,bineq,Aeq,beq,xlb,xub)

nx = size(K,2);
nb = length(b);
clear prob;
[r, res] = mosekopt('symbcon echo(0)');

% Specify the non-conic part of the problem.
prob.c   = [1;sparse(nx,1)];
prob.a   = sparse(0,nx+1);
prob.a   = [prob.a;sparse([zeros(size(Aineq,1),1),Aineq;zeros(size(Aeq,1),1),Aeq;zeros(size(Aeq,1),1),-Aeq])];
prob.buc = sparse([bineq;beq;-beq]);
prob.blx = [-Inf;sparse(xlb)];
prob.bux = [Inf;sparse(xub)];

% Specify the cones as affine conic constraints.
prob.accs = [res.symbcon.MSK_DOMAIN_QUADRATIC_CONE nb+1];
prob.f = sparse([1,zeros(1,nx);zeros(nb,1),K]);
prob.g =  [0;-b];

% Optimize the problem. 
[~,res]=mosekopt('minimize echo(0)',prob);

% Display the primal solution.
sol = res.sol.itr.xx;
x_opt = sol(2:(nx+1));
cost_opt = sol(1)^2;

end