%% Usecase 5 : Optimization with ranked solutions
%   author        : Mirko Pasquini
%
%   In this usecase-example we show how to run a "ranked" optimization,
%   meaning that to each local solution of the optimization we will
%   associate a score, depending on the predicted objective function value,
%   distance from a given dataset and sensitivity with respect to medium
%   variability. Please refer to the "Beginner's guide" for an explanation
%   of these different components.
%
%   In the variable lambda_ranking the user can provide different weights
%   associated to the importance assigned to each one of these score
%   components.
%
%   >>>> IMPORTANT : Execute while being in the usecase-examples folder, 
%   due to relative pathnames. <<<<

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
% Select the type of objective function considered by the optimization. The
% options are 'harvest', 'rate-max' and 'rate-min', respectively maximizing
% the harvest of a product of interest, maximizing its uptake-secretion
% rate or minimizing it.
objective_type = 'harvest';

% Variables for evaluating the distance of the obtained optimized solutions
% from a given dataset of media (in this case contained in medium_temp)
distance_evaluation_vars.norm_type = 2;
distance_evaluation_vars.normalization_flag = 1;
distance_evaluation_vars.media_data = cin_data(index_data.decision_metabolite_indices,:);

% Variables for evaluating the sensitivity of the obtained optimized
% solutions (in particular the sensitivity of associated extracellular
% metabolites concentration and rates with respect to perturbation in the
% medium)
sensitivity_evaluation_vars.met_names =   {'S1','S2','S3','S4','P1','P3','X'};
sensitivity_evaluation_vars.uncertainty_range = [2;2;2;2];
sensitivity_evaluation_vars.N_realizations = 50;

% Weights for the different components that make up the overall score of
% the obtained local solutions of the optimization. In particular the score
% will be affected by objective function values, distance of the local
% solutions from a given media dataset as well as a sensitivity score
% (which in this case is obtained as the sample based standard deviation of
% the predicted harvest)
lambda_ranking.objective = 1;
lambda_ranking.distance = 0.1;
lambda_ranking.sensitivity = 1;

verbose = 1; % verbose = 1 : shows output; verbose = 0 : hides output;


%% Optimization of the product of interest harvest

[optimization_results_harvest, optimization_results_ranked_original] = model_based_optimization_with_ranking(model, ...
                                                                        initial_conditions_data, ...
                                                                        index_data, ...
                                                                        constraints, ...
                                                                        objective_type,...
                                                                        distance_evaluation_vars, ...
                                                                        sensitivity_evaluation_vars, ...
                                                                        lambda_ranking, ...
                                                                        verbose);
%% Re-ranking when importance of objective function is increased

lambda_ranking.objective = 2;
optimization_results_harvest_reranked = get_optimization_results_ranking(optimization_results_harvest, lambda_ranking);

%% Maximization of the product of interest secretion rate

objective_type = 'rate-max';
[optimization_results_rate, optimization_results_ranked_rate] = model_based_optimization_with_ranking(model, ...
                                                                    initial_conditions_data, ...
                                                                    index_data, ...
                                                                    constraints, ...
                                                                    objective_type,...
                                                                    distance_evaluation_vars, ...
                                                                    sensitivity_evaluation_vars, ...
                                                                    lambda_ranking, ...
                                                                    verbose);
                                                                                             

