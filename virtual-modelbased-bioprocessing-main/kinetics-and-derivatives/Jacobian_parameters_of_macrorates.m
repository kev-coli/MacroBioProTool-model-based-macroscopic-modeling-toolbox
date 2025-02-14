%% Jacobian of macrorates with respect to the kinetic parameters
%   author        : Mirko Pasquini
%
%   The function returns the sensitivity of the macrorates with respect to
%   the kinetic parameters, in a matrix form. 
%   The vector of parameters, with respect to which the sensitivity is 
%   evaluated, is obtained by transformation of the parameters matrix into 
%   a vector form (through the function 
%   construct_vector_from_parameters_matrix.m). 
%   Please refer to the "Beginner's guide" for further info.
%
%   JthW = JACOBIAN_PARAMETERS_OF_MACRORATES(theta_matrix,c)
%
%   @inputs:
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       c : vector of extracellular metabolite concentrations (> 0). These 
%           corresponds to coupled metabolites (Please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite")
%
%   @outputs:
%       JthW : sensitivity matrix of the macroreactions with respect to the 
%           kinetic parameters. Due to its structure, the i-th row of the 
%           matrix will have all zero columns except for the indices related
%           to the parameters of the i-th macroreaction.
%           Please refer to the "Beginner's guide" for further information on
%           how the kinetic parameters vector is structured.


function JthW = Jacobian_parameters_of_macrorates(theta_matrix,c)
    m = size(theta_matrix,1); % number of macroreactions
    l = size(theta_matrix,2); % number of parameters per macroreaction
    JthW = zeros(m,m*l);
    HDC_matrix = BuildHDCMatrix(theta_matrix,c); % the HDC matrix is a matrix in which the elements
                                                 % contain evaluation of double components
                                                 % functions. Used for quick access and
                                                 % efficiency
    for i = 1 : m
        nabla_theta_wi = build_ith_macrorate_sensitivity_vector(HDC_matrix(i,:),theta_matrix(i,:),c);
        JthW(i,(i-1)*l+1:i*l) = nabla_theta_wi;
    end
end



