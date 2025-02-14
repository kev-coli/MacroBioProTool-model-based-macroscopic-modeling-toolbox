%% Usecase 4 : Sensitivity analysis and medium distance
%   author        : Mirko Pasquini
%
%   This usecase provides an example of two functionality which are usually
%   connected to the nominal optimization:
%   1. How to perform a sensitivity analysis (both numerical and locally 
%       analytical) to determine how perturbations of the medium
%       composition would affect the prediction of some quantities of
%       interest (i.e. concentrations, rates and harvest);
%   2. How to evaluate the distance of a medium from a dataset of media, to
%       compare the exploratory capacity of the optimization solutions
%       (e.g. an optimized medium that is "more distant" from the data
%       would provide richer information for subsequent modelling
%       purposes).
%
%   Execute while being in the usecase-examples folder, due to relative
%   pathnames.

clear all
close all
clc

%%   Step 0: load the model to use for the optimization and define the file
%%   information for the dataset.

supplementary_file = '../models/K-Net/supplementary.xls';
parameters_file = '../models/K-Net/parameters.xls';

F = 1;          
Xv = 1;  

model = load_model(supplementary_file, parameters_file,F,Xv);

%%  Step 1: define the optimization structures
%   These structure will contain all the information useful in the
%   optimization, such as: constraints, index of the metabolite whose
%   harvest rate we are optimizating, indexes of the decision metabolites,
%   initial conditions etc.

% index_data : contains all the information on indexes relative to the 
%   vector of metabolites, useful for the optimization. This order is the
%   one of the macroreaction stoichiometric matrix Amac (found in the 
%   supplementary.xls file of the model).
%   In this example the metabolites in the concentrations vector are
%   ordered as:
%   ['S1','S2','S3','S4','P1','P3','X']

index_data.decision_metabolite_indices = colvec([1:4]); % i.e. ['S1','S2','S3','S4']
index_data.coupled_met_indices = colvec([1:4]); %i.e. all metabolites except biomass and mAb
index_data.rate_index = 6; % we are interested in maximizing the harvest rate of mAb             
index_data.growth_index = 7; % growth rate index is the rate of biomass

% Lower and upper bounds on the medium composition concentrations. The
% first column refers to the lower bound, the second column refers to 
% the upper bound.
constraints.medium_bounds = [1,   12;
                             1,  12;
                             1, 12;
                             1, 12]; 

% Lower and upper bounds on the metabolites concentrations. Concentrations
% are bounded to be non-negative, and only the concentrations of S2 (index
% 2) and S3 (index 3) are upper-bounded, since they are assumed to be
% toxic byproducts.
constraints.concentration_bounds = [zeros(4,1),[100;
                                                 3;
                                                 4;
                                                 100]]; 


% Solubility constraints. The constraint is a_solub'*u <= b_solub, which in
% this case reads as [S2]_in + [S3]_in <= 18.
constraints.a_solub = [0;1;1;0];
constraints.b_solub = 18;

%% Step 2 : load the data and define initial conditions of optimization
% Load the data from past experiments.
cext_data = readtable("../data/K-Net data/fake_data_example_cext.xlsx");
cext_data = cext_data{2:2:end,3:end}'; % takes last day of each condition
qext_data = readtable("../data/K-Net data/fake_data_example_qext.xlsx");
qext_data = qext_data{2:2:end,3:end}';
cin_data = readtable("../data/K-Net data/fake_data_example_cin.xlsx");
cin_data = cin_data{2:2:end,3:end}'; % takes last day of each condition

% Lower and upper bounds on the metabolites rates. Constraints on the rates
% are not strictly necessary, but allow for the optimizer to remain "close"
% to the data that generated the model. In fact we enforce the rates (with
% the exception of the rate of P3) to stay within 50% of the rates
% measurements in the data.
constraints.rate_bounds = define_rates_bounds_percentage_on_data(qext_data,0.5);
constraints.rate_bounds(6,2) = 100;

% Select all the conditions in the dataset as initial conditions for the 
% optimization.
initial_conditions_data.medium = cin_data(:,1:end);
initial_conditions_data.cext = cext_data(:,1:end);
initial_conditions_data.qext = qext_data(:,1:end);
%% Run optimization from initial conditions from the data
% Execute the optimization algorithm from the initial conditions specified.
% Please refer to the documentation of
% nominal_optimization_given_metabolic_model for a description of
% the output structure containing the optimization results.
optim_results_init_from_data = nominal_optimization_given_metabolic_model(model, initial_conditions_data, index_data, constraints,'iter','harvest');

%% Evaluate sensitivity of quantities of interest
% Here we evaluate the sensitivity of optimal medium 4 i.e. how the
% steady-state concentrations, rates and harvest are affected by a
% variation of such medium.

%% 1. Numerical Sensitivity
%  The numerical sensitivity is evaluated by generating N_realizations (in
%  this case 50) media realizations around the nominal one, and for each of 
%  these media a virtual experiment is executed to generate predictions for
%  steady-state concentrations and rates. After the sensitivity analysis we
%  are able to check mean, std, min and max of steady-state concentrations,
%  rates and harvest with respect to this set of realizations.
met_names =   {'S1','S2','S3','S4','P1','P3','X'};
% indexes of metabolites that are not decision variables of the optimization
not_decision_metabolites = setdiff(index_data.coupled_met_indices, index_data.decision_metabolite_indices); 
nominal_medium = optim_results_init_from_data.media_optimized(:,1);
nominal_concentration = optim_results_init_from_data.cext_optimized(:,1);
uncertainty_range = [2;2;2;2];
N_realizations = 50;
numerical_sensitivity = numerical_sensitivity_medium_variability( model, ...
                                                                met_names, ...
                                                                nominal_medium, ...
                                                                nominal_concentration, ...
                                                                uncertainty_range, ...
                                                                N_realizations, ...
                                                                cin_data(not_decision_metabolites,1), ...
                                                                index_data,...
                                                                1);

%% 2. Local Sensitivity
%  The local sensitivity is evaluated as the gradient of the considered
%  quantities (i.e. steady-state concentrations, rates, and harvest) with
%  respect to small variations in the input. This is a local property, and
%  gives an information of how any individual medium component affect the
%  considered quantity.
local_sensitivity = local_sensitivity_medium_variability(model, ...
                                                         met_names, ...
                                                         nominal_medium, ...
                                                         nominal_concentration,...
                                                         index_data);


%% Evaluate distance of optimized media from dataset
% Here we will evaluate how distant are the optimized media from the
% available data. The idea is to implement solutions that are "different"
% from what has already been implemented.

n_optimized_conditions = size(optim_results_init_from_data.media_optimized,2);
distances = zeros(n_optimized_conditions, 1);
nearest_index = zeros(n_optimized_conditions, 1);

for n = 1 : n_optimized_conditions
    [distances(n), nearest_index(n)] = evaluate_distance_medium_from_dataset(optim_results_init_from_data.media_optimized(:,n),...
                                                                            cin_data(index_data.decision_metabolite_indices,:), ...
                                                                            1, ...
                                                                            1, ...
                                                                            constraints.medium_bounds);
end
