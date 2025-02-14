%% Reorder data based on stoichiometry matrix order
%   author  :   Mirko Pasquini 
%
%   The code takes an unordered data matrix, with the metabolites labels 
%   associated to such data matrix, and reorder the data so that the
%   metabolite order of the macroreaction stoichiometric matrix (specified
%   by the metabolites labels of such matrix) is respected.
%
%   data_reordered = REORDER_DATA_BASED_ON_STOICHIOMETRY(data_unordered, stoichiometry_metabolites_labels, data_metabolites_labels)
%
%   @inputs:
%       data_unordered : data matrix where each column corresponds to a 
%           different tested condition and each row corresponds to a 
%           particular extracellular metabolite. The order of metabolites 
%           here could be different than the order of metabolites in the
%           macroreaction stoichiometric matrix.
%       stoichiometry_metabolites_labels : metabolites labels from the 
%           stoichiometric matrix. The labels are in the form of a cell 
%           array of strings 
%           e.g. txt = {'Glcext'},{'Glnext'},{'Serext'},...
%       data_metabolites_labels  metabolites labels from the 
%           data matrix. The labels are in the form of a cell 
%           array of strings 
%           e.g. txt = {'Glcext'},{'Glnext'},{'Serext'},...
%
%   @outputs:
%       data_reordered : data matrix reordered based on the stiochiometric
%           matrix order

function data_reordered = reorder_data_based_on_stoichiometry(data_unordered, stoichiometry_metabolites_labels, data_metabolites_labels)
    [~, loc_order] = ismember(stoichiometry_metabolites_labels, data_metabolites_labels);
    data_reordered = data_unordered;
    ind_data = 1;
    for k = 1 : length(loc_order)
        if loc_order(k) <= size(data_unordered,1)
            data_reordered(ind_data,:) = data_unordered(loc_order(k),:);
            ind_data = ind_data + 1;
        end
    end
end