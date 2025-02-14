%% Allows to extract the vector of coupled metabolites indexes
%   author  :   Mirko Pasquini
%
%   The function allows to extract the vector of indexes of the coupled
%   metabolites in the vector of all metabolites involved in the
%   stoichiometric matrix (See Beginner's guide for a definition of coupled
%   metabolites).
%
%   coupled_met_indices = EXTRACT_COUPLED_MET_INDICES(parameters_file, stoichiometry_file)
%
%   @inputs: 
%       supplementary_file  : xls file containing information on the
%           macroreaction stoichiometric matrix (on its first sheet called
%           Sheet 1. For the particular structure of the file see an example 
%           file in the "models" folder. Note: please follow this structure 
%           exactly for compatibility.
%       parameters_file : xls file containing the matrix of parameters
%           for the macroreaction kinetics. For the particular structure of 
%           the file see an example file in the "models" folder. Note: please
%           follow this structure exactly for compatibility.
%
%   @outputs:
%       coupled_met_indices :  vector of indexes for the coupled metabolites.
%           The order refers to the order of metabolites in the Amac matrix. 

function coupled_met_indices = extract_coupled_met_indices(parameters_file, stoichiometry_file)
    stoichiometry_labels = readtable(stoichiometry_file,'Sheet','Sheet1');
    stoichiometry_labels = stoichiometry_labels{:,1};
    kinetic_labels = readtable(parameters_file);
    kinetic_labels = kinetic_labels.Properties.VariableNames(2:2:end);
    coupled_met_indices = zeros(length(kinetic_labels),1);
    for i = 1 : length(kinetic_labels)
        for j = 1 : length(stoichiometry_labels)
            if strcmp(stoichiometry_labels{j},kinetic_labels{i})
                coupled_met_indices(i) = j;
            end
        end
    end
end