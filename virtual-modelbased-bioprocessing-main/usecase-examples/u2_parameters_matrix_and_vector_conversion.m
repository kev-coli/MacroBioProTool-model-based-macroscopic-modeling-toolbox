%% Usecase 2 : Convert a parameters matrix into a vector and viceversa
%   author  :   Mirko Pasquini
%
%   This usecase provides an example of how to convert between a (metabolism kinetics)
%   parameters matrix to a parameters vector. We refer to function
%   construct_matrix_from_parameters_vector and
%   construct_vector_from_parameters_matrix for details on how such
%   matrices and vectors are structured.
%
%   The utility of this is not immediate to the end user, but the two
%   structures for the kinetic parameters are used separately in different 
%   functions.
%
%   In this usecase: we load a model (and corresponding paramters matrix),
%   we convert the matrix into vector form and the vector back into matrix
%   form. We show that these operations are one the inverse of the other.
%
%   Execute while being in the usecase-examples folder, due to relative
%   pathnames.

%% Step 1. Load a model

supplementary_file = '../models/K-Net/supplementary.xls';
parameters_file = '../models/K-Net/parameters.xls';

reordered_parameters_pathname = '../models/K-Net/parameters.xls';
% generate_reordered_parameters_file(supplementary_file, parameters_file, reordered_parameters_pathname);
F = 1;          % == Units missing ==
Xv = 60*10.81;  % == Units missing ==

model = load_model(supplementary_file, reordered_parameters_pathname,F,Xv);


%% Step 2. Convert parameters matrix into its vector form

matrix_form_initial = model.theta_matrix;
vector_form = construct_vector_from_parameters_matrix(matrix_form_initial);

%% Step 3. Re-convert parameters vector into its matrix form

matrix_form_final = construct_matrix_from_parameters_vector(vector_form, size(matrix_form_initial,1));

%% Step 4. Shows that the two operations are bijective

infinity_norm_of_matrix_norm_difference = norm(matrix_form_initial - matrix_form_final,Inf); % It is 0 if (and only if) the two matrix forms are equal
