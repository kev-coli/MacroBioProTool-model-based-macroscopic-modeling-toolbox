% This code will run the four steps detailed in MAIN_MODELING.m
% It is launched from MAIN_MODELING.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% STEP 1 : File treatment %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All the options for Step 1 chosen by the user are saved in the structure options_data
all_directories = {directory_file_cell_specific_rate_data_qext;
                   directory_file_concentration_data_cext;
                   directory_file_metabolic_network};

options_data.smoothing = smoothing;
options_data.coeff_smoothing = coeff_smoothing;
options_data.average_based_normalization = average_based_normalization;
options_data.user_defined_normalization_matrix = user_defined_normalization_matrix;
options_data.average_data_per_condition = average_data_per_condition; 
options_data.prediction_media = prediction_media;

% Function for the treatment of the data (STEP 1)
[data_qext,data_cext,data_stoich] = Step_1_file_treatment(all_directories,options_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% STEP 2 : EFM computation %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All the options for Step 2 are saved in the structure options_EFM
options_EFM.use_mosek = use_mosek; 
options_EFM.tolerance_EFM = tolerance_EFM; 
options_EFM.computation_method_EFM = computation_method_EFM;
options_EFM.reduction_of_EFM = reduction_of_EFM;
options_EFM.factor_error_reduction = factor_error_reduction;

% In order to plot some intermediate results, the same specifications as STEP 4 are considered
plot_specification.number_of_plots_per_columns = number_of_plots_per_columns;
plot_specification.number_of_plots_per_rows = number_of_plots_per_rows;

% Function for the EFM computation (STEP 2)
[EFMs,Macro_rates,q_CG_train,q_CG_predict] = Step_2_EFM_computation(data_qext,data_stoich,options_EFM,options_data);
legend_label = ["Model after Step 2"];  % Label used for plots in Step 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STEP 3 : Kinetic identification %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All the options for Step 3 are saved in the structure options_kinetic_identification
options_kinetic_identification.average_data_per_condition = average_data_per_condition;
options_kinetic_identification.mets_not_included_in_Monod_kinetics = mets_not_included_in_Monod_kinetics;
options_kinetic_identification.regularization = regularization;
options_kinetic_identification.lambda_grid = lambda_grid;
options_kinetic_identification.number_multistart_points = number_multistart_points;
options_kinetic_identification.number_max_of_iterations_EM = number_max_of_iterations_EM;
options_kinetic_identification.burn_in_period_EM = burn_in_period_EM;
options_kinetic_identification.number_samples_per_iteration_EM = number_samples_per_iteration_EM;
options_kinetic_identification.number_trials_sampling_EM = number_trials_sampling_EM;
options_kinetic_identification.perturbation_variance_EM = perturbation_variance_EM;
options_kinetic_identification.threshold_variation_neutral_effect = threshold_variation_neutral_effect;
options_kinetic_identification.threshold_fit_activation_or_inhibition = threshold_fit_activation_or_inhibition;

% Function for the Monod kinetic identification (STEP 3)
[Kinetic_parameters,q_model_train,q_model_predict,w_model_train] = Step_3_Monod_kinetic_identification_of_macroscopic_rates(Macro_rates,data_cext,data_qext,data_stoich,EFMs,options_kinetic_identification,options_data);
q_all_models_train = cat(3,q_CG_train,q_model_train);
q_all_models_predict = cat(3,q_CG_predict,q_model_predict);
legend_label = [legend_label;"Model after Step 3"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% STEP 4 : Display results %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_specification.legend_label = legend_label;
[set_of_macro_reations_ext,set_of_macro_reations_meas,all_reactions_involved_in_each_EFM,errors_CG_and_model] = Step_4_treat_results(q_all_models_train,q_all_models_predict,data_qext,data_cext,plot_specification,EFMs,data_stoich,Kinetic_parameters,options_kinetic_identification,options_EFM,options_data,name_save_file);




