%% Test unit script
%
%   Data ---> Model ----> Predictions ----> Plots
%               |-------> Optimization
%
%   Change the parameters for optimization, simulation and modelling, as well
%   as the datafiles, for testing.

clear all
close all
clc

%% Step 1 : run the test unit main function
% ------------------------ START MODIFY FROM HERE -------------------------
% 1.1. load model -- to remove as model should be evaluated in test unit

supplementary_file = '../models/T2024-1-Net2-VN4/supplementary.xls';
parameters_file = '../models/T2024-1-Net2-VN4/parameters.xls';
reordered_parameters_pathname = '../models/T2024-1-Net2-VN4/reordered_parameters.xls';
F = 1;          % == Units missing ==
Xv = 60*10.81;  % == Units missing ==

model = load_model(supplementary_file, reordered_parameters_pathname,F,Xv);

% ------------------------------------------------------------------------
% 1.2. load modelling options
modelling_options = [];
operation_parameters = [];

% ------------------------------------------------------------------------
% 1.3. load datafiles
datafiles.qext_data_file = '../data/qext_avg_datafile.xlsx';
datafiles.qext_data_sheet = 'qext_average_data';
datafiles.cext_data_file = '../data/cext_avg_datafile.xlsx';
datafiles.cext_data_sheet = 'cext_average_data';
datafiles.cin_data_file = '../data/cin_avg_datafile.xlsx';
datafiles.cin_data_sheet = 'cin_average_data';

% ------------------------------------------------------------------------
% 1.4 load optimization options
optimization_options.index_data.decision_metabolite_indices = [1;3;4;12]; % i.e. [Glc, Ser, Asn, Leu]
optimization_options.index_data.coupled_met_indices = [1:23]; %i.e. all metabolites except biomass and mAb
optimization_options.index_data.rate_index = 25; % we are interested in maximizing the harvest rate of mAb             
optimization_options.index_data.growth_index = 24; % growth rate index is the rate of biomass
optimization_options.constraints.medium_bounds = [55,   65;
                             5.1,  20;
                             3.78, 13;
                             3.99, 12]; 
optimization_options.constraints.concentration_bounds = [zeros(23,1),[ones(14,1)*10000;
                                                 20;
                                                 ones(3,1)*10000;
                                                 6;
                                                 ones(4,1)*10000]]; 
[qext_data_avg,~,~] = xlsread(datafiles.qext_data_file,datafiles.qext_data_sheet);
optimization_options.constraints.rate_bounds = define_rates_bounds_percentage_on_data(qext_data_avg',0.35);
optimization_options.constraints.rate_bounds(24,2) = 0.45; %  0.45 for High Growth rate constraint; 0.25 for Low Growth rate constraint
optimization_options.constraints.rate_bounds(25,2) = 100; % mAb rate is not upper bounded.
optimization_options.constraints.a_solub = [0;1;1;0];
optimization_options.constraints.b_solub = 22;
% ------------------------------------------------------------------------
% Here initial conditions for the optimization are taken from actual
% data. Change this part to change initial conditions for the optimization
% (if unchanged than the files with the data should be excel files)
% ------------------------------------------------------------------------
[qext_temp,~,~] = xlsread(datafiles.qext_data_file,datafiles.qext_data_sheet);
[cext_temp,~,~] = xlsread(datafiles.cext_data_file,datafiles.cext_data_sheet);
[medium_temp,~,~] = xlsread(datafiles.cin_data_file,datafiles.cin_data_sheet);
flags.transpose_flag = true;
if flags.transpose_flag 
    qext_temp = qext_temp';
    cext_temp = cext_temp';
    medium_temp = medium_temp';
end
optimization_options.initial_conditions.medium = medium_temp(:,1:15);
optimization_options.initial_conditions.cext = cext_temp(:,1:15);
optimization_options.initial_conditions.qext = qext_temp(:,1:15);
% ------------------------------------------------------------------------
optimization_options.objective_type = 'harvest';
optimization_options.distance_evaluation_vars.norm_type = 2;
optimization_options.distance_evaluation_vars.normalization_flag = 1;
optimization_options.distance_evaluation_vars.media_data = medium_temp(optimization_options.index_data.decision_metabolite_indices,:);
optimization_options.sensitivity_evaluation_vars.met_names =   {'Glc','Gln','Ser','Asn','Asp','Arg','Tyr','Thr','Lys','Val','Ile',...
                'Leu','Phe','Met','Lac','Ala','Gly','Glu','NH4','Cys','His','Pro',...
                'Trp','Biomass','mAb'};
optimization_options.sensitivity_evaluation_vars.uncertainty_range = [3;1.5;1.5;1.5];
optimization_options.sensitivity_evaluation_vars.N_realizations = 100;
optimization_options.lambda_ranking.objective = 1;
optimization_options.lambda_ranking.distance = 0.1;
optimization_options.lambda_ranking.sensitivity = 1;

optimization_options.verbose = 1;

% ------------------------------------------------------------------------
% 1.5 run main function

[model, ...
    optimization_results, ...
    optimization_results_ranked, ...
    model_concentrations_predictions_mbe,...
    model_rates_predictions_mbe,...
    model_rates_predictions_direct,...
    qext_data,...
    cext_data,...
    cin_data] = whole_pipeline_testing(datafiles,...
                                                    modelling_options,...
                                                    operation_parameters,...
                                                    optimization_options,...
                                                    model);
% ------------------------ STOP MODIFY FROM HERE -------------------------

%% Step 2 : draw plots
clc
close all
% ------------------------------------------------------------------------
% 2.1 Draw predictions against actual data

n_metabolites_total = length(optimization_options.sensitivity_evaluation_vars.met_names);
n_metabolites_concentrations = length(optimization_options.index_data.coupled_met_indices);
n_metabolites_per_figure = 9;
n_cols = ceil(sqrt(n_metabolites_per_figure));
n_rows = ceil(n_metabolites_per_figure/n_cols);
n_conditions = size(qext_data,2);
figure_opened = 0;
data_conditions_names = {'N1','N2','N3','N4','N5','N6','N8','N9','N10','N11','N12','N13','N14',...
                            'N15','N17','T1','T2','T3','T4','T5','T6','T7','T8','T9'};

% ------------------------------------------------------------------------
% MBE based predictions: predictions obtained by solving the
% mass-balance equation (input is the medium composition)

% mbe concentrations predictions
for n = 1 : n_metabolites_concentrations
    figure(ceil(n/n_metabolites_per_figure));
    sgtitle('MBE Concentrations predictions vs real data')
    n_sp = mod(n,n_metabolites_per_figure);
    if n_sp == 0 
        n_sp = n_metabolites_per_figure;
    end
    subplot(n_rows, n_cols, n_sp);
    hold on
    grid on
    plot(1:n_conditions,cext_data(n,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,model_concentrations_predictions_mbe(n,:),'-ob','LineWidth',2.0);
    xticks([1:n_conditions]);
    xticklabels(data_conditions_names);
    xlabel('Condition');
    ylabel(strcat('[',optimization_options.sensitivity_evaluation_vars.met_names{n},'] ([mmol/L])'));
end

figure_opened = ceil(n/n_metabolites_per_figure);

% mbe rates predictions
for n = 1 : n_metabolites_total
    figure(figure_opened + ceil(n/n_metabolites_per_figure));
    sgtitle('MBE Rates predictions vs real data')
    n_sp = mod(n,n_metabolites_per_figure);
    if n_sp == 0 
        n_sp = n_metabolites_per_figure;
    end
    subplot(n_rows, n_cols, n_sp);
    hold on
    grid on
    plot(1:n_conditions,qext_data(n,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,model_rates_predictions_mbe(n,:),'-ob','LineWidth',2.0);
    xlabel('Condition');
    xticks([1:n_conditions]);
    xticklabels(data_conditions_names);
    ylabel(strcat('q_',optimization_options.sensitivity_evaluation_vars.met_names{n}, '([1/day])'));
end

figure_opened = figure_opened + ceil(n/n_metabolites_per_figure);
% ------------------------------------------------------------------------
% direct rates predictions: predictions based directly on the kinetic model
% (inputs are the measured concentrations in the data)
for n = 1 : n_metabolites_total
    figure(figure_opened + ceil(n/n_metabolites_per_figure));
    sgtitle('Direct Rates predictions vs real data')
    n_sp = mod(n,n_metabolites_per_figure);
    if n_sp == 0 
        n_sp = n_metabolites_per_figure;
    end
    subplot(n_rows, n_cols, n_sp);
    hold on
    grid on
    plot(1:n_conditions,qext_data(n,:),'-xk','LineWidth',2.0);
    plot(1:n_conditions,model_rates_predictions_direct(n,:),'-ob','LineWidth',2.0);
    xlabel('Condition');
    xticks([1:n_conditions]);
    xticklabels(data_conditions_names);
    ylabel(strcat('q_',optimization_options.sensitivity_evaluation_vars.met_names{n}, '([1/day])'));
end

figure_opened = figure_opened + ceil(n/n_metabolites_per_figure);

% ------------------------------------------------------------------------
% 2.2 Draw plots for optimized solution with sensitivity

% concentrations
optimization_results_ranked = get_optimization_results_ranking(optimization_results, optimization_options.lambda_ranking);

n_optimized_conditions = size(optimization_results_ranked.media_optimized,2);
optimized_conditions_names = {'M1','M2','M3','M4','M5',...
                              'M6','M8','M9','M10','M11',...
                              'M12','M13','M14','M15','M17'};


sensitivities_concentrations = zeros(n_metabolites_concentrations,n_optimized_conditions);
sensitivities_rates = zeros(n_metabolites_total,n_optimized_conditions);
sensitivities_harvest_rates = zeros(n_metabolites_total,n_optimized_conditions); % not all harvest rates are of interest, but only the one of the POI

for n = 1 : n_metabolites_total
    for k = 1 : n_optimized_conditions
        if n <= (n_metabolites_concentrations)
            sensitivities_concentrations(n,k) = optimization_results_ranked.sensitivity_structures{k}{n}.concentrations.std;
        end
        sensitivities_rates(n,k) = optimization_results_ranked.sensitivity_structures{k}{n}.rates.std;
        sensitivities_harvest_rates(n,k) = optimization_results_ranked.sensitivity_structures{k}{n}.harvest.std;
    end
end


for n = 1 : n_metabolites_concentrations
    figure(figure_opened + ceil(n/n_metabolites_per_figure));
    sgtitle('Optimized solutions predictions -- concentrations')
    n_sp = mod(n,n_metabolites_per_figure);
    if n_sp == 0 
        n_sp = n_metabolites_per_figure;
    end
    subplot(n_rows, n_cols, n_sp);
    hold on
    grid on
    errorbar(1:n_optimized_conditions,...
                optimization_results.cext_optimized(n,:),...
                sensitivities_concentrations(n,:),'-ob','LineWidth',2.0);
    xlabel('Condition');
    xticks([1:n_optimized_conditions]);
    xticklabels(optimized_conditions_names);
    ylabel(strcat('[',optimization_options.sensitivity_evaluation_vars.met_names{n},'] ([mmol/L])'));
end

figure_opened = figure_opened + ceil(n/n_metabolites_per_figure);

% rates

for n = 1 : n_metabolites_total
    figure(figure_opened + ceil(n/n_metabolites_per_figure));
    sgtitle('Optimized solutions predictions -- rates')
    n_sp = mod(n,n_metabolites_per_figure);
    if n_sp == 0 
        n_sp = n_metabolites_per_figure;
    end
    subplot(n_rows, n_cols, n_sp);
    hold on
    grid on
    errorbar(1:n_optimized_conditions,...
            optimization_results.qext_optimized(n,:),...
            sensitivities_rates(n,:),...
            '-ob','LineWidth',2.0);
    xlabel('Condition');
    xticks([1:n_optimized_conditions]);
    xticklabels(optimized_conditions_names);
    ylabel(strcat('q_',optimization_options.sensitivity_evaluation_vars.met_names{n}, '([1/day])'));
end

figure_opened = figure_opened + ceil(n/n_metabolites_per_figure);

% harvest
figure(figure_opened + 1);
title('Optimized solutions predictions -- harvest rate')
hold on
grid on
errorbar(1:n_optimized_conditions,...
        optimization_results.objective_function_value_optimized,...
        sensitivities_harvest_rates(optimization_options.index_data.rate_index,:),...
        '-ob','LineWidth',2.0);
xlabel('Condition');
xticks([1:n_optimized_conditions]);
xticklabels(optimized_conditions_names);
ylabel('Harvest rate of POI');




