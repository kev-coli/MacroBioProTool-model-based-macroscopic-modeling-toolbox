%% Solve the nominal optimization problem (explicit formulation)
%   author : Mirko Pasquini 
%
%   The nominal optimization problem that we want to solve is
%
%   maximize{u,c}   (F-q_growth)*q_k               (Harvest rate)          
%   subject to      g1: u - u_upper <= 0            (Medium bounds)
%                   g2: -u + u_lower <= 0
%                   g3: c - c_upper <= 0            (Concentration bounds)
%                   g4: c + c_upper <= 0            
%                   g5: q - q_upper <= 0            (Rate bounds)
%                   g6: -q + q_lower <= 0
%                   g7: a'u - b <= 0                (Solubility constraint)
%                   h : Xv*q - F*c + F[u;u_hat] = 0 (Mass Balance Eq.)
%
%   where q are the extracellular metabolite rates (with q_growth being the 
%   growth rate and q_k being the rate of the metabolite whose harvest we 
%   aim to maximize), c are the extracellular metabolite concentrations, 
%   while u are the decision variable of the optimization, i.e. the 
%   concentrations of the extracellular metabolites that we can change in 
%   the medium (the ones that are fixed are denoted by u_hat).
%
%   Remark: 
%   we refer to this formulation as explicit, because we have both
%   medium (u) and concentrations (c) as decision variables, which are linked
%   together by the last equality constraint (the mass-balance equation).
%
%     [z,lambda,fval,exitflag,first_ord] = ...
%       SOLVE_NOMINAL_MEDIUM_OPTIMIZATION_EXPLICIT(model, u0, c0, u_bounds,...
%       c_bounds, q_bounds, a_solub, b_solub, coupled_met_indices,...
%       display_option, u_hat, rate_index, growth_index)
%
%   @inputs
%       model : this is a structure containing all the useful information on 
%           the kinetic model. Please refer to the "Beginner's guide" for 
%           further info.
%       u0 : vector of initial medium composition, relative to the decision
%           metabolites, used as initial condition for the optimization.
%           Please refer to the "Beginner's guide" for a definition of 
%           decision metabolite.
%       c0 : vector of initial coupled metabolite concentration, used as 
%           initial condition for the optimization. Please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite. 
%       u_bounds : bounds for the medium composition (only for metabolites 
%           that are modified in the medium). Matrix with two columns, the 
%           first column is the vector of lower bounds, while the second
%           column is the vector of upper bounds;
%       c_bounds : bounds for the extracellular metabolite concentrations. 
%           Matrix with two columns, the first column is the vector of lower 
%           bounds, while the second column is the vector of upper bounds;
%       q_bounds : bounds for the extracellular metabolite rates. Matrix 
%           with two columns, the first column is the vector of lower 
%           bounds, while the second column is the vector of upper bounds;
%       a_solub : vector 'a' of the solubility linear constraint
%       b_solub : scalar 'b' of the solubility linear constraint
%       coupled_met_indices : vector of indexes for the coupled metabolites. 
%           The order refers to the order of metabolites in the model.Amac 
%           matrix. Please refer to the "Beginner's guide" for a definition
%           of coupled metabolite.
%       display_option : see optimoptions for fmincon.
%       u_hat : vector containing the components of the medium that are 
%           fixed, i.e. the ones that are not relative to decision variables 
%           of the optimization problem
%       rate_index : index of the rate of interest in the uptake/secretion 
%           rates vector
%       growth_index : index of the growth rate in the uptake/secretion     
%           rates vector
%       max_iterations: specify the maximum number of iterations for the
%           optimization. If this is not specified the default number will
%           be 1e3.
%       objective_type: specify the type of objective to optimize. The
%           available options are:
%           'harvest' : maximize the predicted harvest of a certain 
%               metabolite of interest (whose index is rate_index)
%           'rate-min' : minimize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%           'rate-max' : maximize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%
%   @outputs
%       z : concatenation of the medium vector and the vector of 
%           extracellular metabolite concentrations, obtained as a result
%           of the optimization procedure.
%       lambda : structure containing the Lagrangian multipliers of the
%           optimization
%           lambda.g1 : multipliers associated to the upper bound
%               constraints on the decision variable z
%           lambda.g2 : multipliers associated to the lower bound
%               constraints on the decision variable z
%           lambda.g3 : multipliers associated to the upper bound 
%               constraints on the extracellular metabolite rates q
%           lambda.g4 : multipliers associated to the lower bound 
%               constraints on the extracellular metabolite rates q
%           lambda.g5 : multipliers associated to the solubility constraint
%           lambda.h : multipliers associated to the mass-balance equation
%               equality constraints
%       fval : objective function value obtained as a result of the
%           optimization. It is negative in the case of maximization (e.g.
%           harvest rate)
%       exitflag : exit flag of the optimization procedure. Please refer to
%           fmincon documentation for a description of the exit flags. 
%           As a rule of thumb if the exit flag is positive then the
%           optimization can be considered succesful.
%       first_ord : First-order optimality criterion obtained at the end of
%           the optimization procedure. Please refer to fmincon documentation
%           for a description of the first-order optimality criteria. As a
%           rule of thumb if this value is close to 0 we can consider
%           that the optimization procedure converged to a feasible local 
%           optimizer (we consider any value below 1e-4 to be acceptable).
%
% ========================================================================
%
%       NOTE: the optimization assumes that the components of the 
%       medium vector that are optimized, are in the first positions 
%       of the medium vector (i.e. the metabolite order is taken such 
%       that the medium vector is [u;u_hat]). Please consider calling 
%       the function reorder_model_for_optimization on the model 
%       structure, before using it in the nominal optimization. 
% 
%       NOTE: In general this function should not be called directly, 
%       instead the interface nominal_optimization_given_metabolic_model 
%       should be considered, where such re-ordering is not required.
%
% ========================================================================

function [z,lambda,fval,exitflag,first_ord] = solve_nominal_medium_optimization_explicit(model,u0,c0,u_bounds,c_bounds,q_bounds,a_solub,b_solub,coupled_met_indices,display_option,u_hat,rate_index,growth_index,max_iterations, objective_type)
    fmincon_options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Display',display_option,'CheckGradients',false,'MaxIterations', max_iterations);
    z0 = [u0;c0];
    [z_ip,fval_ip,exitflag_ip,output_ip,lambda_ip] = fmincon(@(z)objective_medium_optimization(z,model,length(u0),rate_index,growth_index, objective_type),z0,[a_solub',zeros(1,length(c0))],b_solub,[],[],[u_bounds(:,1);c_bounds(:,1)],[u_bounds(:,2);c_bounds(:,2)],@(z)nonlinear_constraints_medium_optimization(z,model,q_bounds,length(u0),coupled_met_indices,u_hat),fmincon_options);
    lambda_ip = define_lambda_structure(lambda_ip,q_bounds);
    % optimization with the interior-point is usually better at converging
    % to an optimum, but the sqp method returns more reliable Lagrangian
    % multipliers (lambda)
    fmincon_options = optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'Display',display_option,'CheckGradients',false,'MaxIterations', max_iterations);
    z0 = z_ip;
    [z_sqp,fval_sqp,exitflag_sqp,output_sqp,lambda_sqp] = fmincon(@(z)objective_medium_optimization(z,model,length(u0),rate_index,growth_index, objective_type),z0,[a_solub',zeros(1,length(c0))],b_solub,[],[],[u_bounds(:,1);c_bounds(:,1)],[u_bounds(:,2);c_bounds(:,2)],@(z)nonlinear_constraints_medium_optimization(z,model,q_bounds,length(u0),coupled_met_indices,u_hat),fmincon_options);
    lambda_sqp = define_lambda_structure(lambda_sqp,q_bounds);
    if exitflag_sqp >= 0
        z = z_sqp;
        lambda = lambda_sqp;
        exitflag = exitflag_sqp;
        first_ord = output_sqp.firstorderopt;
        fval = -fval_sqp; % it is a maximization problem, which is cast in minimization of
                          % minus the objective function.
    else
        z = z_ip;
        lambda = lambda_ip;
        exitflag = exitflag_ip;
        first_ord = output_ip.firstorderopt;
        fval = -fval_ip; % it is a maximization problem, which is cast in minimization of
                          % minus the objective function.
    end
end

function lambda  = define_lambda_structure(lambda_temp,q_bounds)
    nq = size(q_bounds,1);
    lambda.g1 = lambda_temp.upper;
    lambda.g2 = lambda_temp.lower;
    lambda.g3 = lambda_temp.ineqnonlin(1:nq);
    lambda.g4 = lambda_temp.ineqnonlin(nq+1:end);
    lambda.g5 = lambda_temp.ineqlin;
    lambda.h = lambda_temp.eqnonlin;
end

function [fobj,gradf] = objective_medium_optimization(z,model,nu,rate_index,growth_index,objective_type)
    Amac = model.Amac;
    Xv = model.Xv;
    F = model.F;
    theta = model.theta_matrix;
    u = z(1:nu);
    c = z(nu+1:end);
    q = ComputeQExt(c,theta,Amac);
    interested_macrorates = (model.Amac(growth_index,:)~=0) + (model.Amac(rate_index,:)~=0);
    JcQ = Amac*MacroKineticsJacobian(theta,c,interested_macrorates);
    if strcmp(objective_type, 'harvest') % maximize harvest of product of interest
        fobj = -(F-q(growth_index))*q(rate_index);
        gradf = -[zeros(nu,1);
                 -JcQ(growth_index,:)'*q(rate_index)+(F-q(growth_index))*JcQ(rate_index,:)'];
    elseif strcmp(objective_type, 'rate-min') % minimize rate of interest
        fobj = q(rate_index);
        gradf = [zeros(nu,1);
                 JcQ(rate_index,:)'];
    else % maximize rate of interest
        fobj = -q(rate_index);
        gradf = -[zeros(nu,1);
                 JcQ(rate_index,:)'];
    end
end

function [cin,ceq,dc,dceq] = nonlinear_constraints_medium_optimization(z,model,q_bounds,nu,coupled_met_indices,u_fixed)
    Amac = model.Amac;
    Xv = model.Xv;
    F = model.F;
    theta = model.theta_matrix;
    nc = length(z)-nu;
    u = z(1:nu);
    c = z(nu+1:end);
    q = ComputeQExt(c,theta,Amac);
    nq = length(q);
    JcQ = Amac*MacroKineticsJacobian(theta,c);
    cin = [q - q_bounds(:,2);
         -q + q_bounds(:,1)];
    dc = [zeros(nu,2*nq);
          [JcQ', -JcQ']];
    F1 = F*[zeros(nu,nc-nu);eye(nc-nu)];
    F2 = F*[eye(nu);zeros(nc-nu,nu)];
    ceq = Xv*q(coupled_met_indices)-F*c+F2*u+F1*u_fixed;
    dceq = [F2';
            Xv*JcQ(coupled_met_indices,:)' - F*eye(nc)];
end