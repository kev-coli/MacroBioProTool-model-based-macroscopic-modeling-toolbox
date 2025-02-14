%% Preprocessing of parameters file to re-order metabolites
%   author  :   Mirko Pasquini 
%
%   This function is used to preprocess the parameters file, so that the
%   metabolites are in the same order as the stoichiometric matrix in the
%   supplementary file.
%
%   @inputs: 
%       supplementary_file : xls file containing information on the
%           macroreaction stoichiometric matrix on its first sheet. For the
%           particular structure of the file see an example file in the
%           "models" folder. Note: please follow this structure exactly for
%           compatibility.
%       parameters_file : xls file containing the matrix of parameters
%           for the macroreaction kinetics. For the particular structure of
%           the file see an example file in the "models" folder. Note:
%           please follow this structure exactly for compatibility.
%       new_parameters_file_pathname : complete pathname string for the
%           reordered parameters xls file (i.e. it contains the filename)
%
%   @outputs:
%       reordering_map : it is a vector which remap the order of the 
%           metabolites in the parameters matrix to the order of the 
%           metabolites in the stoichiometric matrix (i.e. if metabolite X 
%           is in position i in the parameters file, while it is in position
%           j of the stoichiometric matrix file, then reordering_map(i)=j) 

function reordering_map = generate_reordered_parameters_file(supplementary_file, parameters_file, new_parameters_file_pathname)
    [~,txt,~] = xlsread(supplementary_file,'Sheet1');
    [or_para,txt_par,~] = xlsread(parameters_file);
    [reordered_theta_matrix, reordering_map] = reorder_parameters_based_on_stoichiometry(txt,txt_par(1,2:2:end),or_para);
    % generate reordered labels for the parameters
    reordered_txt_par = txt_par;
    for i = 1 : length(reordering_map)
        reordered_txt_par(1,2*reordering_map(i)) = txt_par(1,2*i);
        reordered_txt_par(2,2*reordering_map(i)) = txt_par(2,2*i);
        reordered_txt_par(1,2*reordering_map(i)+1) = txt_par(1,2*i+1);
        reordered_txt_par(1,2*reordering_map(i)+1) = txt_par(1,2*i+1);        
    end
    % generate a cell matrix to convert in an excel file
    new_file_cell_array = cell(size(reordered_theta_matrix,1)+size(txt_par,1),size(txt_par,2));
    new_file_cell_array(1:2,:) = reordered_txt_par;
    for i = 1 : size(reordered_theta_matrix,1)
        for j = 1 : size(reordered_theta_matrix,2)
            new_file_cell_array{i+2,j} = reordered_theta_matrix(i,j);
        end
    end
    % generate the new excel file
    xlswrite(new_parameters_file_pathname,new_file_cell_array);
end