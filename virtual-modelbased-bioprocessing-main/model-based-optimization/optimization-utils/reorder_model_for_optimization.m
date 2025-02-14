%% Re-order model components to match the optimization metabolite order 
%   author  :   Mirko Pasquini
%
%   This function is used to reorder the rows of the macroreaction
%   stoichiometric matrix, and the columns of the kinetic parameters
%   matrix, so that the decision metabolites are placed at the beginning
%   of the metabolites vector (hence the consequent re-arrangement of rows
%   and columns for Amac and theta_matrix respectively). Please refer to the
%   "Beginner's guide" for a definition of decision metabolite.
%
%   NOTE: THE END USER SHOULD NOT BE CONCERNED BY THIS FUNCTION IN GENERAL.
%   TO PERFORM NOMINAL OPTIMIZATION PLEASE USE THE APPROPRIATE INTERFACES,
%   WHICH DO NOT REQUIRE A RE-ORDERING OF THE METABOLITES VECTOR.
%
%   model_reordered = REORDER_MODEL_FOR_OPTIMIZATION(model,decision_metabolites_indices,coupled_met_indices)
%   
%   @inputs:
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide"
%           for further info.
%       decision_metabolites_indexes : vector of indexes of the metabolites 
%           in the medium vector that has to be optimized (i.e. these are 
%           the main decision variables of the optimization problem). Please
%           refer to the "Beginner's guide" for a definition of decision 
%           metabolite.
%       coupled_met_indexes : vector of indexes for the coupled metabolites. 
%           The order refers to the order of metabolites in the original 
%           Amac matrix. Please refer to the "Beginner's guide" for a 
%           definition of coupled metabolite.
%
%   @outputs:
%       model_reordered : this is a structure containing all the useful 
%           information on the re-ordered kinetic model. Please refer to the
%           "Beginner's guide" for further info. In this model the first 
%           rows of Amac and the first columns of theta_matrix corresponds
%           to decision metabolites.

function model_reordered = reorder_model_for_optimization(model,decision_metabolites_indexes,coupled_met_indexes)
    Amac = model.Amac;
    theta = model.theta_matrix;
    not_decision_metabolites_indices = setdiff(1:size(Amac,1),decision_metabolites_indexes); % all metabolites for which rates are expressed
    Amac_dec_var_rows = Amac(decision_metabolites_indexes,:); % rows of Amac relative to decision metabolites
    Amac_not_dec_var_rows = Amac(not_decision_metabolites_indices,:); % rows of Amac relative to non-decision metabolites
    not_decision_metabolites_indices = setdiff(coupled_met_indexes,decision_metabolites_indexes); % only metabolites that appear in the rates expressions.

    % columns of theta_matrix relative to decision metabolites
    theta_dec_var_columns = [];
    for i = 1 : length(decision_metabolites_indexes)
        theta_dec_var_columns = [theta_dec_var_columns,theta(:,decision_metabolites_indexes(i)*2:decision_metabolites_indexes(i)*2+1)];
    end
    
    % columns of the_matrix relative to non-decision metabolites
    theta_not_dec_var_columns = [];
    for i = 1 : length(not_decision_metabolites_indices)
        theta_not_dec_var_columns = [theta_not_dec_var_columns,theta(:,not_decision_metabolites_indices(i)*2:not_decision_metabolites_indices(i)*2+1)];
    end
    Amac_reordered = [Amac_dec_var_rows;
                      Amac_not_dec_var_rows];

    theta_reordered = [theta(:,1), theta_dec_var_columns, theta_not_dec_var_columns];
    model_reordered = model;
    model_reordered.Amac = Amac_reordered;
    model_reordered.theta_matrix = theta_reordered;
end
