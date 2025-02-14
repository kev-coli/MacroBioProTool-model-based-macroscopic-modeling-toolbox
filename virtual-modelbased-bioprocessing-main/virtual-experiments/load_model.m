%% Load model structure from supplementary and parameters files
%   author        : Mirko Pasquini 
%
%   The function takes as inputs the pathname of the supplementary and
%   parameters files, as well as the flow rate F and the viable cells at
%   steady state Xv, and load the model structure that will be used by many
%   other functions in the library.
%
%   NOTE: In this version we assume that supplementary_file and
%   parameters_file (as well as any data file) have the metabolites in the
%   same order. If this is not the case these files should be
%   pre-processed.
%
%   model = LOAD_MODEL(supplementary_file, parameters_file, F, Xv)
%
%   @inputs: 
%       supplementary_file  : xls or mat file containing information on the
%           macroreaction stoichiometric matrix (on its first sheet if an xls
%           file). For the particular structure of the file see an example 
%           file in the "models" folder. Note: please follow this structure 
%           exactly for compatibility.
%       parameters_file : xls or mat file containing the matrix of parameters
%           for the macroreaction kinetics. For the particular structure of 
%           the file see an example file in the "models" folder. Note: please
%           follow this structure exactly for compatibility.
%       F : flow rate of perfusion 
%       Xv : viable cells in the bioreactor at steady-state 
%
%   @outputs:
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide" 
%           for further info.

function model = load_model(supplementary_file, parameters_file, F, Xv)
    [~,~,supplementary_file_extension] = fileparts(supplementary_file);
    if (strcmp(supplementary_file_extension,".xlsx") || strcmp(supplementary_file_extension,".xls"))
        Amac = readtable(supplementary_file,'Sheet','Sheet1');
        Amac = Amac{:,2:end};
    elseif (strcmp(supplementary_file_extension,".mat"))
        Amac = load(supplementary_file);
    else
        disp("Unknown extension of supplementary file. Only '.xlsx', '.xls' and '.mat' are allowed.");
    end

    [~,~,parameters_file_extension] = fileparts(parameters_file);
    if (strcmp(parameters_file_extension,".xlsx") || strcmp(parameters_file_extension,".xls"))
        theta_matrix = readtable(parameters_file);
        theta_matrix = theta_matrix{:,:};
    elseif (strcmp(parameters_file_extension,".mat"))
        theta_matrix = load(parameters_file);
    else
        disp("Unknown extension of supplementary file. Only '.xlsx', '.xls' and '.mat' are allowed.");
    end

    model.theta_matrix = theta_matrix;
    model.Amac  = Amac;
    model.F = F;
    model.Xv = Xv;
end