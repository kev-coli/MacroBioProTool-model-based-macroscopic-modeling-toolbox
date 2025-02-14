%% Jacobian of extracellular metabolites rates with respect to the kinetic parameters
%   author        : Mirko Pasquini 
%
%   The function returns the sensitivity of the extracellular metabolites 
%   rates with respect to the kinetic parameters, in a matrix form. 
%   The vector of parameters, with respect to which the sensitivity is 
%   evaluated, is obtained by transformation of the parameters matrix into 
%   a vector form (through the function 
%   construct_vector_from_parameters_matrix.m). 
%   Please refer to the "Beginner's guide" for further info.
%
%   JthQ = Jacobian_parameters_of_rates(theta_matrix,c, Amac)
%
%   @inputs:
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       c : vector of extracellular metabolite concentrations (> 0). These 
%           corresponds to coupled metabolites (Please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite")
%       Amac : macroreaction stoichiometric matrix, obtained as Aext*E, 
%           where Aext is the extracellular stoichiometric matrix and E 
%           is the matrix of EFMs.
%
%   @outputs:
%       JthQ : sensitivity matrix of the extracellular metabolite rates 
%           with respect to the vector of kinetic parameters. Please refer 
%           to construct_vector_from_parameters_matrix.m for how the total 
%           parameters vector is structured.

function JthQ = Jacobian_parameters_of_rates(theta_matrix,c, Amac)
    JthW = Jacobian_parameters_of_macrorates(theta_matrix,c);
    JthQ = Amac*JthW;
end