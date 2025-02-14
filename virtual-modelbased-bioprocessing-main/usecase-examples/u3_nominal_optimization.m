%% Usecase 3 : Nominal model-based medium optimization for kinetic model
%   author  :   Mirko Pasquini
%
%   This usecase provides an example of how to run the nominal optimization
%   procedure for a given model. We will start the optimization using the 
%   already tested media as initial conditions of the optimization.
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
