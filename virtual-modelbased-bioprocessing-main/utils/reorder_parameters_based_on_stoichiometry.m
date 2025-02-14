%% Reorder the matrix of parameters based on the stoichiometric matrix order
%   author        : Mirko Pasquini 
%
%  Many times the order of metabolites in the parameters matrix and the one 
%  in the macroreaction stoichiometric matrix is  different, and this 
%  causes problems in the simulations. The function takes as inputs the two 
%  vector of metabolite labels (from the parameters file and from the 
%  stoichiometric matrix file) and the original parameters matrix, and 
%  reorder the parameters so that they now match with the stoichiometric 
%  matrix labels.
%
%  Example:
%       Stoichiometric matrix:
%                   w1  w2  w3
%          1: Glc   x   x   x
%          2: Lac   x   x   x
%          3: Asn   x   x   x
%
%       Parameters matrix (original):
%          wmax    Asn     Glc     Asn  
%          wmax    ka  ki  ka  ki  ka  ki
%          W1      x   x   y   y   z   z                 
%          W2      x   x   y   y   z   z  
%          W3      x   x   y   y   z   z
%
%       Parameters matrix (reordered):
%          wmax    Glc     Lac     Asn  
%          wmax    ka  ki  ka  ki  ka  ki
%          W1      y   y   x   x   z   z                 
%          W2      y   y   x   x   z   z  
%          W3      y   y   x   x   z   z
%
%   [theta_matrix_ordered, loc_par] = REORDER_PARAMETERS_BASED_ON_STOICHIOMETRY(txt,txt_par,para)
%
%   @inputs:
%       txt : metabolites labels from the stoichiometric matrix. The labels 
%           are in the form of a cell array of strings 
%           e.g. txt = {'Glcext'},{'Glnext'},{'Serext'},...
%       txt_par : metabolites labels from the parameters matrix. The labels
%           are in the form of a cell array of strings
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide"
%
%   @outputs:
%       theta_matrix_ordered : matrix of kinetic parameters, ordered
%           based on the labels order in the macroreaction stoichiometric 
%           matrix.
%       reordering_map : it is a vector which remap the order of the 
%           metabolites in the parameters matrix to the order of the 
%           metabolites in the stoichiometric matrix (i.e. if metabolite X 
%           is in position i in the parameters file, while it is in position 
%           j of the stoichiometric matrix file, then reordering_map(i)=j) 

function [theta_matrix_ordered, reordering_map] = reorder_parameters_based_on_stoichiometry(txt,txt_par,theta_matrix)
    txt_par = colvec(txt_par);
    i = 1;
    while i <= length(txt)
        if isempty(txt{i})
            txt = txt(i+1:end); % some initial elements of txt can be empty 
                                % and this causes problems in the following 
                                % use of the function. This remove these
                                % elements from txt
        else
            i = i + 1;          % check next element. only here because if
                                % first element is removed now the iterator
                                % points exactly to the new first element
                                % we want to examine
        end
    end

    [~,reordering_map] = ismember(txt_par,txt);
    nr = size(theta_matrix,1);
    n_met = length(txt); %includes mab and biomass generally
    para_ordered_res = reorder_parameters_based_on_stoichiometric_matrix(theta_matrix(:,2:end),reordering_map,n_met,nr);
    theta_matrix_ordered = [theta_matrix(:,1),para_ordered_res];
end

function coupled_metabolites_para_residual_ordered = reorder_parameters_based_on_stoichiometric_matrix(para_residual_unordered,reordering,n_met,nr)
    para_residual_ordered = zeros(nr,n_met*2);
    for i = 1 : length(reordering)
        para_residual_ordered(:,reordering(i)*2-1:reordering(i)*2) = para_residual_unordered(:,i*2-1:i*2);
    end
    % 1. understand which numbers are missing from reordering
    columns_to_remove = setdiff(1:n_met,reordering);
    columns_to_maintain = setdiff(1:n_met,columns_to_remove);
    % 2. remove the corresponding columns from the para_residual_ordered
    % matrix, since those are the columns corresponding to mAb and biomass
    % that do not appear in the parameters matrix
    coupled_metabolites_para_residual_ordered = zeros(nr,length(columns_to_maintain));
    for i = 1 : length(columns_to_maintain)
        k = columns_to_maintain(i);
        coupled_metabolites_para_residual_ordered(:,[2*i-1,2*i]) = para_residual_ordered(:,[2*k-1,2*k]);
    end
end