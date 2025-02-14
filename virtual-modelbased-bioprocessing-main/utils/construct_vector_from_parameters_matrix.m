%% Generate a vector of parameters starting from a matrix
%   author : Mirko Pasquini 
%
%   The function takes a matrix of parameters and output the parameters 
%   in a vector form. This function and the "opposite" one, which is 
%   construct_matrix_from_parameters_vector, allows to put the parameters 
%   in the form that is more useful to the particular considered function.
%
%   theta_vector = CONSTRUCT_VECTOR_FROM_PARAMETERS_MATRIX(theta_matrix)
%
%   @inputs  
%       theta_matrix : matrix of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.
%
%   @outputs
%       theta_vector : vector of macroreaction kinetics parameters. Please 
%           refer to the "Beginner's guide" for further info.


function theta_vector = construct_vector_from_parameters_matrix(theta_matrix)
    m = size(theta_matrix,1);
    q = size(theta_matrix,2);
    theta_vector = zeros(1,m*q);
    for i = 1 : m
        theta_vector((i-1)*q+1:i*q) = theta_matrix(i,:);
    end
    theta_vector = theta_vector';
end