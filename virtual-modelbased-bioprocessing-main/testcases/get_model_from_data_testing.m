%% TESTING : Get model from data testing function
%
%
%   Input : data and operation_parameters (i.e. Xv and F)
%
%   Output : model
%


function model = get_model_from_data_testing(data, modelling_parameters, operation_parameters)
    model = macroscopic_modeling_from_data(data, modelling_parameters, operation_parameters);
end