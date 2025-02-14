%% Build total medium composition vector given fixed and decision metabolites
%   author  :   Mirko Pasquini
%
%   When optimizing, the medium composition vector consists in a fixed part
%   (u_hat) and a variable part (u) which corresponds to the decision
%   variables of the optimization.
%
%   u_total = BUILD_U_TOTAL(u, u_hat, index_data)
%
%   @inputs:
%       u : vector containing the components of the medium relative to
%           decision variables of the optimization problem
%       u_hat : vector containing the components of the medium that are
%           fixed, i.e. the ones that are not relative to decision
%           variables of the optimization problem
%       index_data : a structure containing all information related to 
%           important indexes for the optimization.
%               index_data.decision_metabolite_indices : vector of indexes 
%                   of the metabolites in the medium vector that has to be
%                   optimized (i.e. these are the main decision variables
%                   of the optimization problem).
%               index_data.coupled_met_indices : vector of indexes for the 
%                   coupled metabolites. The order refers to the order of 
%                   metabolites in the Amac matrix. Please refer to the 
%                   "Beginner's guide" for a definition of coupled 
%                   metabolite.
%               index_data.rate_index : index of the rate of interest in the
%                   uptake/secretion rates vector
%               index_data.growth_index : index of the growth rate in the 
%                   uptake/secretion rates vector 
%
%   @outputs:
%       u_total : vector containing both decision and non-decision
%           metabolites, in the right order (following the one in
%           index_data.decision_metabolite_indices).


function u_total = build_u_total(u, u_hat, index_data)
    nu_tot = length([u; u_hat]);
    u_total = zeros(nu_tot,1);
    u_total(index_data.decision_metabolite_indices) = u;
    u_total(setdiff(1:nu_tot,index_data.decision_metabolite_indices)) = u_hat;
end