%% Usecase 1 : Load a model from excel files and run a virtual experiments
%   author  :   Mirko Pasquini
%
%   This usecase provides an example of how a model can be loaded
%   starting from excel files on parameters and stoichiometric matrix.
%   After the model is loaded a virtual experiment is executed.
%
%   As an additional example we simulate the conditions from the data using
%   the loaded model (can be used as a model sanity check).
%
%   IMPORTANT : Execute while being in the usecase-examples folder, due to 
%   relative pathnames.

clc
close all
clear all

%%   Step 0: define some complete pathnames for the various files

supplementary_file = '../models/K-Net/supplementary.xls';
parameters_file = '../models/K-Net/parameters.xls';

%%   Step 1: generate reordered parameters file to match the same metabolite
%           order as the supplementary file. This should be done as the
%           load_model function requires that the metabolites in the
%           supplementary and parameters files are in the same order

reordered_parameters_pathname = '../models/K-Net/parameters.xls';
% generate_reordered_parameters_file(supplementary_file, parameters_file, reordered_parameters_pathname);

%%   Step 2: Load a model using the load_model function. You need to specify
%           the flow rate F and the viable cells Xv

F = 1;          
Xv = 1;  

model = load_model(supplementary_file, reordered_parameters_pathname, F, Xv);


%%   Step 3.a : Run a virtual experiment
%            The model structure contains all the information necessary to
%            run a virtual experiment. To run a virtual experiment we need
%            to specify a vector of medium concentrations (u) and a vector 
%            of initial concentrations (c0). These two vectors have to be
%            of the same size. The order of metabolites is the same as
%            specified in the supplementary excel file. The vector of
%            initial concentrations c0 is just a starting condition for the
%            solver that tries to solve the mass-balance equation.

%   The metabolites in the concentrations vector are ordered as:
%   ['S1','S2','S3','S4','P1','P2','P3','X']
cext_data = readtable("../data/K-Net data/fake_data_example_cext.xlsx");
cext_data = cext_data{2:2:end,3:end}'; % takes last day of each condition
qext_data = readtable("../data/K-Net data/fake_data_example_qext.xlsx");
qext_data = qext_data{2:2:end,3:end}';
cin_data = readtable("../data/K-Net data/fake_data_example_cin.xlsx");
cin_data = cin_data{2:2:end,3:end}'; % takes last day of each condition

c0 = cext_data(:,1);
q = ComputeQExt(cext_data(:,1),model.theta_matrix, model.Amac); %% TEST

coupled_met_indices = 1:4; 
display_option = 'iter-detailed'; % it will show the iterations of the 
                                  % mass-balance equation solver
u = cin_data(:,1);
[c, q, w] = run_virtual_experiment(u, model, c0, coupled_met_indices, display_option);
mbe = model.Xv*q(coupled_met_indices) - model.F*c + model.F*u;

%% Step 3.b : Compare measurements data with model predictions
% We first load the data


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
