%% Generate a matrix of parameters starting from a vector
%   author        : Mirko Pasquini 
%
%   The function takes a vector of parameters and the number of 
%   macroreactions, and output the parameters in a matrix form. This 
%   function and the "opposite" one, which is 
%   construct_vector_from_parameters_matrix, allows to put the parameters 
%   in the form that is more useful to the particular considered function.
%
%   theta_matrix = CONSTRUCT_MATRIX_FROM_PARAMETERS_VECTOR(theta_vector, m)
%
%   @inputs  
%       theta_vector : vector of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%       m : number of macroreactions (this will be the number of rows of
%           theta_matrix).
%
%   @output
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.


function theta_matrix = construct_matrix_from_parameters_vector(theta_vector, m)
    q = length(theta_vector)/m;
    theta_matrix = zeros(m,q);
    for i = 1 : m
        theta_matrix(i,:) = theta_vector((i-1)*q+1:i*q);
    end
end