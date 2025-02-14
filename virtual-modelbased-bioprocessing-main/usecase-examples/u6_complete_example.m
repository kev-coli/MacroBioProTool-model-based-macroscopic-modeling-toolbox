%% Usecase 6 : Complete example from modelling to optimization through simulations
%
%   This usecase provides an example of how the whole library can be used. 
%   First the measurements data are used to generate a kinetic model,
%   then the model is checked with predictions against the actual data, and
%   at the end the model-based medium optimization is run to find the
%   medium the maximize the predicted product harvest.
%
%   IMPORTANT : Execute while being in the usecase-examples folder, due to relative
%   pathnames.

clear all
close all
clc

addpath(genpath('virtual-modelbased-bioprocessing-main')); % This allows to add the pathway to the entire toolbox

%% 1. Model Identification

%   The metabolites in the concentrations vector are ordered as:
%   ['S1','S2','S3','S4','P1','P2','P3','X']
%% WARNING - WARNING
% ALL THE PARAMETERS LINKED TO IDENTIFICATION HAVE TO BE CHANGED 
% IN THE MAIN_MODELING.m SCRIPT IN THE FOLDER macroscopic-modeling

run("../macroscopic-modeling/MAIN_MODELING.m")

load(strcat("../models/",name_save_file,"/",name_save_file,".mat"))
model.Amac = Amac;
model.theta_matrix = cell2mat(Parameters_cell(3:end,:)); 
model.F = 1;
model.Xv = 1;

clc
close all


%% 2. Predictions check
%   Model predictions of the actual data are checked here

cext_data = readtable("../data/K-Net data/fake_data_example_cext.xlsx");
cext_data = cext_data{2:2:end,3:end}'; % takes last day of each condition
qext_data = readtable("../data/K-Net data/fake_data_example_qext.xlsx");
qext_data = qext_data{2:2:end,3:end}';
cin_data = readtable("../data/K-Net data/fake_data_example_cin.xlsx");
cin_data = cin_data{2:2:end,3:end}'; % takes last day of each condition

coupled_met_indices = 1:4; 
n_conditions = size(cext_data,2);
% variables to contain model predictions
cext_model = cext_data;
qext_model = qext_data;

for n = 1 : n_conditions
    fprintf('Simulating condition n. %d /n',n);
    [cext_model(:,n), qext_model(:,n), ~] = run_virtual_experiment(cin_data(:,n),model,cext_data(:,n),coupled_met_indices,'final');
end

figure(1)
    subplot(2,2,1);
    hold on
    grid on
    xlabel('Condition #');
    ylabel('S1');
    plot(1:n_conditions,cext_data(1,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,cext_model(1,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(2,2,2);
    hold on
    grid on
    xlabel('Condition #');
    ylabel('S2');
    plot(1:n_conditions,cext_data(2,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,cext_model(2,:),'-ob','LineWidth',2.0);
    legend('data','model');    
    subplot(2,2,3);
    hold on
    grid on
    xlabel('Condition #');
    ylabel('S3');
    plot(1:n_conditions,cext_data(3,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,cext_model(3,:),'-ob','LineWidth',2.0);
    legend('data','model');    
    subplot(2,2,4);
    hold on
    grid on
    xlabel('Condition #');
    ylabel('S4');
    plot(1:n_conditions,cext_data(4,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,cext_model(4,:),'-ob','LineWidth',2.0);
    legend('data','model');


figure(2)
    subplot(3,3,1)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{S1}');
    plot(1:n_conditions,qext_data(1,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(1,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,2)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{S2}');
    plot(1:n_conditions,qext_data(2,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(2,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,3)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{S3}');
    plot(1:n_conditions,qext_data(3,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(3,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,4)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{S4}');
    plot(1:n_conditions,qext_data(4,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(4,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,5)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{P1}');
    plot(1:n_conditions,qext_data(5,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(5,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,6)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{P3}');
    plot(1:n_conditions,qext_data(6,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(6,:),'-ob','LineWidth',2.0);
    legend('data','model');
    subplot(3,3,7)
    hold on
    grid on
    xlabel('Condition #');
    ylabel('q_{X}');
    plot(1:n_conditions,qext_data(7,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,qext_model(7,:),'-ob','LineWidth',2.0);
    legend('data','model');


%% 3. Model-based optimization

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
