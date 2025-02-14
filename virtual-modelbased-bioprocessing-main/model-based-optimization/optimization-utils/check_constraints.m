%% Check if constraints of the optimization problem are satisfied
%   author  :   Mirko Pasquini 
%
%   This function takes as input a model, the optimization constraints and
%   a particular medium composition and steady-state concentrations, and 
%   return true if the medium and the concentrations satisfy the
%   constraints for the given model, or false otherwise.
%
%   tf = CHECK_CONSTRAINTS(model,u,c,constraints)
%
%   @inputs:
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide" 
%           for further info.
%       u : vector containing the components of the medium relative to 
%           decision variables of the optimization problem
%       c : vector of extracellular metabolite concentrations (> 0). These 
%           corresponds to coupled metabolites (please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite")
%       constraints : a structure containing all information related to the 
%           constraints of the optimization problem (please refer to other
%           functions  such as "model_based_optimization_with_ranking" for
%           further info on this structure)
%
%   @outputs:
%       tf : True if u and c (and the corresponding rates q, given the
%       model) satisfy the optimization constraints specified in
%       constraints. False otherwise.



function tf = check_constraints(model,u,c,constraints)
    q = ComputeQExt(c,model.theta_matrix,model.Amac);
    tf = true;
    if sum(u < constraints.medium_bounds(:,1)) > 0 || sum(u > constraints.medium_bounds(:,2)) > 0
        tf = false; % at least one medium component is outside the feasible bounds
    elseif sum(c < constraints.concentration_bounds(:,1)) > 0 || sum(c > constraints.concentration_bounds(:,2)) > 0
        tf = false; % at least one concentration is outside the feasible bounds
    elseif sum(q < constraints.rate_bounds(:,1)) > 0 || sum(q > constraints.rate_bounds(:,2)) > 0
        tf = false; % at least one uptake-secretion rate is outside the feasible bounds
    elseif constraints.a_solub'*u > constraints.b_solub
        tf = false; % solubility constraint is violated       
    end
end