%% Re-order data rows to match the optimization metabolite order 
%   author  :   Mirko Pasquini 
%
%   This function is used to reorder the rows of the data matrices (i.e. 
%   media compositions and measurements of the corresponding concentrations 
%   and rates), so that the decision metabolites are placed at the beginning
%   of the metabolites vector (hence the consequent re-arrangement of rows 
%   for the data, in which each row correspond to a different metabolite). 
%   Please refer to the "Beginner's guide" for a definition of decision 
%   metabolite.
%
%   NOTE: THE END USER SHOULD NOT BE CONCERNED BY THIS FUNCTION IN GENERAL.
%   TO PERFORM NOMINAL OPTIMIZATION PLEASE USE THE APPROPRIATE INTERFACES,
%   WHICH DO NOT REQUIRE A RE-ORDERING OF THE METABOLITES VECTOR.
%
%   [medium_data, cext_data, qext_data] = 
%       REORDER_DATA_FOR_OPTIMIZATION(medium_temp, cext_temp, qext_temp, decision_metabolite_indices, coupled_met_indices)
%   
%   @inputs:
%       medium_temp : data matrix for media where each column corresponds 
%           to a different tested condition and each row corresponds to a 
%           particular extracellular metabolite. The order of metabolites
%           here is the one corresponding to the original macroreaction
%           stoichiometric matrix (the one to which
%           decision_metabolite_indexes and coupled_metabolite_indexes
%           refer).
%       cext_temp : data matrix for concentrations where each column 
%           corresponds to a different tested condition and each row 
%           corresponds to a particular extracellular metabolite. The order
%           of metabolites here is the one corresponding to the original 
%           macroreaction stoichiometric matrix (the one to which
%           decision_metabolite_indexes and coupled_metabolite_indexes
%           refer).
%       qext_temp : data matrix for uptake-secretion rates where each column 
%           corresponds to a different tested condition and each row 
%           corresponds to a particular extracellular metabolite. The order
%           of metabolites here is the one corresponding to the original 
%           macroreaction stoichiometric matrix (the one to which
%           decision_metabolite_indexes and coupled_metabolite_indexes
%           refer).
%       decision_metabolite_indexes : vector of indexes of the metabolites 
%           in the medium vector that has to be optimized (i.e. these are 
%           the main decision variables of the optimization problem). Please
%           refer to the "Beginner's guide" for a definition of decision 
%           metabolite.
%       coupled_metabolite_indexes : vector of indexes for the coupled 
%           metabolites. The order refers to the order of metabolites in the
%           Amac matrix. Please refer to the "Beginner's guide" for a 
%           definition of coupled metabolite.
%   
%   @outputs:
%       medium_data : data matrix for media where the rows are re-ordered
%           to match the optimization order (i.e. the first rows correspond
%           to decision metabolites)
%       cext_data : data matrix for concentrations where the rows are re-ordered
%           to match the optimization order (i.e. the first rows correspond
%           to decision metabolites)
%       qext_data : data matrix for rates where the rows are re-ordered
%           to match the optimization order (i.e. the first rows correspond
%           to decision metabolites)

function [medium_data, cext_data, qext_data] = reorder_data_for_optimization(medium_temp, cext_temp, qext_temp, decision_metabolite_indexes, coupled_metablite_indexes)
    nq = size(qext_temp,1);
    not_decision_metabolites_indices_for_rates = setdiff(1:nq, decision_metabolite_indexes); % all metabolites for which rates are expressed
    not_decision_metabolites_indices_for_concentrations = setdiff(coupled_metablite_indexes, decision_metabolite_indexes); % only the couple metabolites
    medium_data = [medium_temp(decision_metabolite_indexes,:);
                   medium_temp(not_decision_metabolites_indices_for_concentrations,:)];
    cext_data = [cext_temp(decision_metabolite_indexes,:);
                 cext_temp(not_decision_metabolites_indices_for_concentrations,:)];
    qext_data = [qext_temp(decision_metabolite_indexes,:);
                 qext_temp(not_decision_metabolites_indices_for_rates,:)];
end