%   [set_of_macro_reactions_ext,set_of_macro_reactions_meas,all_reactions_involved_in_each_EFM,errors_CG_and_model] = Step_4_treat_results(q_all_models_train,q_all_models_predict,data_qext,data_cext,plot_specification,EFMs,data_stoich,Kinetic_parameters,options_kinetic_identification,options_EFM,options_data,name_save_file)
%
%   This code treats, plots and saves several modleing results obtained after the three first steps
%
%   @inputs:
%       - q_all_models_train : 3D matrix containing the computed macrosocpic rates from the column generation
%       - q_all_models_predict : 3D matrix containing all the information about the data cext which is used in STEPs 2, 3 and 4
%       - data_cext : structure containing all the information about the data cext which is used in STEPs 2, 3 and 4
%               -> data_cext.ic : same as data_qext.ic
%               -> data_cext.mets_meas_in_cext : cell with all the names of the extracellular metabolites whose concentration is measured
%               -> data_cext.media : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_cext.cext : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_cext.cext_train : matrix of the concentration TRAINING data (averaged, smoothed and/or normalized)
%               -> data_cext.cext_predict : matrix of the concentration PREDICTION data (averaged, smoothed and/or normalized)
%       - data_qext : structure containing all the information about the data qext which is used in STEPs 2, 3 and 4
%               -> data_qext.ia : vector used in order to identify the condition to which each measurement corresponds to (combined with data_qext.ic) 
%               -> data_qext.ic : vector used in order to identify the condition to which each measurement corresponds to (combined with data_qext.ia)  
%               -> data_qext.media : cell vector with the name of each experimental condition identified in the qext qnd cext data files
%               -> data_qext.n_media : number of different experimental conditions
%               -> data_qext.index_training_media : vector gathering all the indices of the conditions in media which are dedicated for model training
%               -> data_qext.index_prediction_media : vector gathering all the indices of the conditions in media which are dedicated for model prediction
%               -> data_qext.N : number of measurements (data)
%               -> data_qext.mets_meas_in_qext : cell with all the names of the extracellular metabolites whose cell specific rate is measured
%               -> data_qext.qext : matrix of the cell specific rate data. Each row corresponds to a given metabolite and the columns correspond to the measurement time instants
%               -> data_qext.qext_train : matrix of the cell specific rate TRAINING data (averaged, smoothed and/or normalized)
%               -> data_qext.qext_predict : matrix of the cell specific rate PREDICTION data (averaged, smoothed and/or normalized)
%               -> data_qext.media_column : cell vector containing the name of the media/condition to which each measurement comes from
%               -> data_qext.days_column : cell vector containing the measurement time instants of each condition
%       - data_stoich : structure containing all the information about the metabolic network which is used in STEPs 2, 3 and 4
%               -> data_stoich.mets : cell vector with all the metabolites involved in the metabolic network   
%               -> data_stoich.A : stoichiometric matrix. The rows correspond to the metabolites and the columns to the reactions 
%               -> data_stoich.n_v : number of reactions. It also corresponds to the number of columns in the stoichiometric matrix A
%               -> data_stoich.is_ext : boolean vector of same dimension as mets. It the j-th entry of is_ext is 1 (resp. 0), then the j-th metabolite in mets is extracellular (resp. intracellular).
%               -> data_stoich.mets_ext : cell vector which gathers the extracellular metabolites 
%               -> data_stoich.mets_int : cell vector which gathers the intracellular metabolites 
%               -> data_stoich.n_mets_int : number of intracellumar metabolites
%               -> data_stoich.rev_bool : boolean vector whose dimension is equal to the number of reactions (n_v). It the j-th entry of rev_bool is 1 (resp. 0), then the j-th reaction is reversible (resp. irreversible).
%               -> data_stoich.ind_rev : vector made up of the indices of the reversible reactions in the flux vector 
%               -> data_stoich.ind_irrev : vector made up of the indices of the irreversible reactions in the flux vector 
%               -> data_stoich.n_v_rev : number of reversible reactions (equal to the dimension of ind_rev)
%               -> data_stoich.n_v_irrev : number of irreversible reactions (equal to the dimension of ind_irrev)
%               -> data_stoich.Dirrev : matrix of dimension n_v_irrev*n_v such that v_irrev = D_irrev*v with v_irrev the subvector of the flux vector v containing only the irreversible fluxes. Useful in order to implement the constraint Dirrev*v >=0. 
%               -> data_stoich.Aext : extracellular stoichiometric matrix. The rows correspond to the EXTRACELLULAR metabolites and the columns to the reactions
%               -> data_stoich.Aint : intracellular stoichiometric matrix. The rows correspond to the  INTRACELLULAR metabolites and the columns to the reactions
%               -> data_stoich.is_meas : boolean vector whose dimension equals the number of extracellular metabolites. If the j-th entry of is_meas is 1 (resp. 0), then the cell specific rate of the j-th extracellular metabolite in mets_ext is measured (resp. not measured)
%               -> data_stoich.Ameas : extracellular stoichiometric matrix with only the measured extrecallular metabolites.
%               -> data_stoich.n_meas : number of extrecallular metabolites whose cell specific rate is measured
%               -> data_stoich.mets_meas_in_stoich : measured metabolites which could be identified in the metabolic network
%       - plot_specification : structure definingthe plot options
%               -> plot_specification.number_of_plots_per_columns : number of plots per columns for displaying the cell specific rate modeling results
%               -> plot_specification.number_of_plots_per_rows : number of plots per rows for displaying the cell specific rate modeling results
%               -> plot_specification.legend_label : cell vector containing the label for the legend of the plots
%       - EFMs : EFM matrix
%       - Kinetic_parameters : structure containing the identified kinetic parameters
%               -> Kinetic_parameters.activation_parameters : matrix containing all the identified activation parameters. The rows correspond to the macro-rates and the column to the concentration (input)
%               -> Kinetic_parameters.inhibition_parameters : matrix containing all the identified inhibition parameters. The rows correspond to the macro-rates and the column to the concentration (input)
%               -> Kinetic_parameters.max_reaction_rate : vector containing all the identified maximal rate constants of each macroscopic rate.
%               -> Kinetic_parameters.lambda_l : optimal regularization parameter selected by the cross-validation. If no regularization was desired, then lambda_l is set to 0.
%       - options_kinetic_identification : structure containing all the options for the Monod kinetic identification
%               -> options_kinetic_identification.average_data_per_condition : boolean. If it is equal to 1 (resp. 0), all the data of the same ocndition are averaged (resp. not averaged). 
%               -> options_kinetic_identification.mets_not_included_in_Monod_kinetics : cell vector made up of all the metabolites which should not be included in the Monod kinetics model
%               -> options_kinetic_identification.regularization : boolean. To add regularization, regularization should be set to 1. 
%               -> options_kinetic_identification.lambda_grid : vector of non-negative scalars. It contains the different values for the regularization parameters which will be tried for model selection (ONLY IF regularization is set to 1)
%               -> options_kinetic_identification.number_multistart_points : number of initial parameter optimizations realized during Step 3.a with different initial conditions. The more initial identifications, the greater the chances for the computation of the global optimum
%               -> options_kinetic_identification.number_max_of_iterations_EM : number of maximal iterations that will be done for the hyperparameter tuning during the Expectation Maximization algorithm
%               -> options_kinetic_identification.burn_in_period_EM : burn-in period of the Expectation Maximization. It corresponds to a warm-up period of the algorithm in order to guarantee proper hyperparemeter tuning. It should be taked large enough.
%               -> options_kinetic_identification.number_samples_per_iteration_EM : number of parameter samples used for the tuning of the hyperparameters
%               -> options_kinetic_identification.number_trials_sampling_EM : number of trials allowed for the sampling with the Metropolis-Hastings algorithm
%               -> options_kinetic_identification.perturbation_variance_EM : non-negative scalar. It corresponds to add additive term added to the estimate of the standard deviation of the log-normal prior distribution. It allows a larger exploration of the kinetic parameter space for the sampling.
%               -> options_kinetic_identification.threshold_variation_neutral_effect : threshold for reduction of an identified modulation function into a neutral effect. It should be non-negative and kept close to 0%
%                                                                                the closer to 0%, the less the chances to possibly reduce some identfiied functions into a neutral effect.
%               -> options_kinetic_identification.threshold_fit_activation_or_inhibition : threshold for reduction of an identified modulation function into an activation or inhibition effect. It should be non-negative and kept close to 100% 
%                                                                                          the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%       - options_EFMs : structure containing all the options for the EFM computation :
%               -> options_EFM.use_mosek : boolean variable. If equal to 1, the column generation algorithm will use MOSEK. IT MUST BE INSTALLED AND THE PATHWAY INDICATED IN MAIN_MODELIN.
%               -> options_EFM.tolerance_EFM : threshold under which the modeling error for the EFM computation is judged negligible enough : all the most relevant EFMs are therefore computed
%               -> options_EFM.computation_method_EFM : character chain. It indicates the chosen method used for the computation of the EFMs. It can also be equal to 'global', 'sequential' and 'union'.
%               -> options_EFM.reduction_of_EFM : boolean. It specifies if the EFM should be reduced after their computation with the column generation (if equal to 1, reduction happens. If taken equal to 0, no reduction is done).
%               -> options_EFM.factor_error_reduction : least-squares error degration allowed for EFM reduction.
%       - options_data : a structure with three entries
%               -> options_data.smoothing : boolean. 0 if the user does not want to smooth the data cext and qext, 1 if he does
%               -> options_data.coeff_smoothing : span for the smoothing as used in the Matlab fonction smooth.
%               -> options_data.normalization : boolean. 0 if the user does not want to normalize the data qext, 1 if he does
%               -> options_data.normalization_matrix : normalization matrix defined by the user. If it is empty, the code will consider a diagonal matrix with the average of the data in the diagonal elements
%               -> options_data.average_data_per_condition : boolean. If it is equal to 1 (resp. 0), all the data of the same ocndition are averaged (resp. not averaged). 
%               -> options_data.prediction_media : cell vector containing the name of the conditions used for prediction
%       - name_save_file : the name of the folder and the .mat and .xls files where the model information will be saved
%
%   @outputs:
%       - set_of_macro_reations_ext : cell vector containing the macroscopic reactions with the extracellular metabolites (measured or not)
%       - set_of_macro_reations_meas : cell vector containing the macroscopic reactions with the MEASURED extracellular metabolites only
%       - all_reactions_involved_in_each_EFM : cell vector of diimension equal to the number of macroscopic reactions. 
%                                              Each cell contains the reactions of the metabolite network which are active in a given macroscopic reaction/EFM
%       - errors_CG_and_model : structure containing absolute and relative error between the training and prediction data and the model without kinetics and with kinetics
%               -> errors_CG_and_model.relative_error_training_qext_CG_model : vector made up of all the relative errors of each cell specific rate between the training data and the model without kinetic         
%               -> errors_CG_and_model.relative_error_training_qext_kinetic_model : vector made up of all the relative errors of each cell specific rate between the training data and the model with kinetic 
%               -> errors_CG_and_model.relative_error_prediction_qext_CG_model : vector made up of all the relative errors of each cell specific rate between the prediction data and the model without kinetic         
%               -> errors_CG_and_model.relative_error_prediction_qext_kinetic_model : vector made up of all the relative errors of each cell specific rate between the prediction data and the model with kinetic 
%               -> errors_CG_and_model.absolute_error_training_qext_CG_model : vector made up of all the absolute errors of each cell specific rate between the training data and the model without kinetic         
%               -> errors_CG_and_model.absolute_error_training_qext_kinetic_model : vector made up of all the absolute errors of each cell specific rate between the training data and the model with kinetic 
%               -> errors_CG_and_model.absolute_error_prediction_qext_CG_model : vector made up of all the absolute errors of each cell specific rate between the prediction data and the model without kinetic         
%               -> errors_CG_and_model.absolute_error_prediction_qext_kinetic_model : vector made up of all the absolute errors of each cell specific rate between the prediction data and the model with kinetic 

function [set_of_macro_reactions_ext,set_of_macro_reactions_meas,all_reactions_for_each_EFM,errors_CG_and_model,Amac] = Step_4_treat_results(q_all_models_train,q_all_models_predict,data_qext,data_cext,plot_specification,EFMs,data_stoich,Kinetic_parameters,options_kinetic_identification,options_EFM,options_data,name_save_file)

  clc
  fprintf("Step 4: Computation of the relative and aboslute error for the measured cell specific rates qext...")
  pause(2) 
  
  % Load the cell specific rate training and prediction data
  qext_train = data_qext.qext_train';
  qext_predict = data_qext.qext_predict';
  nb_data_per_training_media = data_qext.nb_data_per_training_media;
  nb_data_per_prediction_media = data_qext.nb_data_per_prediction_media;

  % Load the properties of the stoichiometric matrix
  mets_meas_in_qext = data_qext.mets_meas_in_qext;
  mets_ext = data_stoich.mets_ext;
  Aext = data_stoich.Aext;
  A = data_stoich.A;
  mets = data_stoich.mets;

  % Computation of the training and prediction RELATIVE error for each cell specific rate for the column generation model (without kinetics) and the kinetics model 
  relative_error_training_qext_CG_model = sum(((qext_train - q_all_models_train(:,:,1)).^2)./(mean(qext_train,2)),2)/sum(nb_data_per_training_media);
  relative_error_training_qext_kinetic_model = sum(((qext_train - q_all_models_train(:,:,2)).^2)./(mean(qext_train,2)),2)/sum(nb_data_per_training_media);
  errors_CG_and_model.relative_error_training_qext_CG_model = relative_error_training_qext_CG_model;
  errors_CG_and_model.relative_error_training_qext_kinetic_model = relative_error_training_qext_kinetic_model;
  if(size(qext_predict,2) > 0)
    relative_error_prediction_qext_CG_model = sum(((qext_predict - q_all_models_predict(:,:,1)).^2)./(mean(qext_predict,2)),2)/sum(nb_data_per_prediction_media);
    relative_error_prediction_qext_kinetic_model = sum(((qext_predict - q_all_models_predict(:,:,2)).^2)./(mean(qext_predict,2)),2)/sum(nb_data_per_prediction_media);
  else
    relative_error_prediction_qext_CG_model = [];  
    relative_error_prediction_qext_kinetic_model = [];  
  end
  errors_CG_and_model.relative_error_prediction_qext_CG_model = relative_error_prediction_qext_CG_model;
  errors_CG_and_model.relative_error_prediction_qext_kinetic_model = relative_error_prediction_qext_kinetic_model;
 
  % Computation of the training and prediction ABSOLUTE error for each cell specific rate for the column generation model (without kinetics) and the kinetics model 
  absolute_error_training_qext_CG_model = sum(((qext_train - q_all_models_train(:,:,1)).^2),2)/sum(nb_data_per_training_media);
  absolute_error_training_qext_kinetic_model = sum(((qext_train - q_all_models_train(:,:,2)).^2),2)/sum(nb_data_per_training_media);  
  errors_CG_and_model.absolute_error_training_qext_CG_model = absolute_error_training_qext_CG_model;
  errors_CG_and_model.absolute_error_training_qext_kinetic_model = absolute_error_training_qext_kinetic_model;
  if(size(qext_predict,2) > 0)
    absolute_error_prediction_qext_CG_model = sum(((qext_predict - q_all_models_predict(:,:,1)).^2),2)/sum(nb_data_per_prediction_media);
    absolute_error_prediction_qext_kinetic_model = sum(((qext_predict - q_all_models_predict(:,:,2)).^2),2)/sum(nb_data_per_prediction_media);
  else
    absolute_error_prediction_qext_CG_model = [];  
    absolute_error_prediction_qext_kinetic_model = [];
  end
  errors_CG_and_model.absolute_error_prediction_qext_CG_model = absolute_error_prediction_qext_CG_model;
  errors_CG_and_model.absolute_error_prediction_qext_kinetic_model = absolute_error_prediction_qext_kinetic_model;
  fprintf("Done.\n") 
  pause(1)

  % Display the table with the computed relative errors
  fprintf("\n")
  fprintf("Relative error on the cell specific rates (qext).\n \n")
  if(size(qext_predict,2) > 0)
    Table_relative_error = table(mets_meas_in_qext',relative_error_training_qext_CG_model,relative_error_prediction_qext_CG_model,relative_error_training_qext_kinetic_model,relative_error_prediction_qext_kinetic_model);
    Table_relative_error.Properties.VariableNames = ["Metabolite","Relative training error (CG)","Relative prediction error (CG)","Relative training error (kinetic)","Relative prediction error (kinetic)"];
  else
    Table_relative_error = table(mets_meas_in_qext',relative_error_training_qext_CG_model,relative_error_training_qext_kinetic_model);
    Table_relative_error.Properties.VariableNames = ["Metabolite","Relative training error (column generation)","Relative training error (kinetic model)"];  
  end  
  disp(Table_relative_error)
  pause(2) 
  
  % Display the table with the computed absolute errors
  fprintf("\n")
  fprintf("Absolute error on the cell specific rates (qext).\n \n")
  if(size(qext_predict,2) > 0)
    Table_absolute_error = table(mets_meas_in_qext',absolute_error_training_qext_CG_model,absolute_error_prediction_qext_CG_model,absolute_error_training_qext_kinetic_model,absolute_error_prediction_qext_kinetic_model);
    Table_absolute_error.Properties.VariableNames = ["Metabolite","Absolute training error (CG)","Absolute prediction error (CG)","Absolute training error (kinetic)","Absolute prediction error (kinetic)"];
  else
    Table_absolute_error = table(mets_meas_in_qext',absolute_error_training_qext_CG_model,absolute_error_training_qext_kinetic_model);
    Table_absolute_error.Properties.VariableNames = ["Metabolite","Absolute training error (CG)","Absolute training error (kinetic)"];  
  end
  disp(Table_absolute_error)
  pause(2) 
  
  fprintf("Computing the macroscopic reactions and pathways...") 
  pause(2)
  fprintf("Done.\n") 
  Aext_mac = Aext*EFMs; % Extracellular macroscopic stoichiometric matrix
  Aext_mac(abs(Aext_mac) < 1e-8) = 0; % Remove negligible term
  
  Ameas = data_stoich.Ameas;
  Amac = Ameas*EFMs;
  Amac(abs(Amac) < 1e-8) = 0;
  
  % Computation of the macroscopic reactions with all the extracellular metabolites. We also compute the reactions which are involved in each EFM
  [set_of_macro_reactions_ext,all_reactions_for_each_EFM] = compute_macro_reactions(Aext_mac,mets_ext,A,mets,EFMs); 
  
  % Computation of the macroscopic reactions with only the MEASURED extracellular metabolites.
  [set_of_macro_reactions_meas,~] = compute_macro_reactions(Amac,mets_meas_in_qext,A,mets,EFMs);
  pause(2) 

  fprintf("Plotting the results for the cell specific rates...\n")
  pause(2) 
  legend_label = plot_specification.legend_label;

  number_of_plots_per_rows = plot_specification.number_of_plots_per_rows;
  number_of_plots_per_columns = plot_specification.number_of_plots_per_columns;
  
  % Creating the x-axis label with the name of the media for the plot
  media = data_qext.media;  
  index_training_media = data_qext.index_training_media;
  index_prediction_media = data_qext.index_prediction_media;
  x_axis_start_ticks = [1,cumsum(nb_data_per_training_media) + 1]; last_tick_training = x_axis_start_ticks(end);
  x_axis_start_ticks = [x_axis_start_ticks,cumsum(nb_data_per_prediction_media) + last_tick_training];
  x_axis_start_ticks(end) = [];
  x_axis_label = [media(index_training_media);media(index_prediction_media)];
 
  % Matlab function which plots the data and the output of the column generation model and kinetic model for comparison. This is done for both training and prediction data
  plot_qext_data_vs_model(qext_train,q_all_models_train,qext_predict,q_all_models_predict,data_qext,x_axis_start_ticks,x_axis_label,number_of_plots_per_rows,number_of_plots_per_columns,legend_label)

  %% Save the necessary variables for model-based optimization
  directory_for_saving = strcat("../models/",name_save_file);
  mkdir(directory_for_saving); % create a new folder with the name of the model 

  fprintf("Saving the files required for model-based optimization in the folder %s...\n",directory_for_saving)
  pause(2)     

  % Macroscopic stoichiometric matrix with only the MEASURED extracellular metabolites
  Ameas = data_stoich.Ameas;
  Amac = Ameas*EFMs;
  Amac_cell = cell(size(Amac,1),size(Amac,2)+1);
  Amac_cell(:,1) = data_stoich.mets_meas_in_stoich;
  Amac(abs(Amac) < 1e-8) = 0;
  Amac_cell(:,2:end) = num2cell(Amac);
  delete(strcat(directory_for_saving,"/supplementary.xls"))
  writecell(Amac_cell,strcat(directory_for_saving,"/supplementary.xls"))
  
  % Macroscopic stoichiometric matrix with all the extracellular metabolites 
  Aext_mac_cell = cell(size(Aext_mac,1),size(Aext_mac,2)+1);
  Aext_mac_cell(:,1) = data_stoich.mets_ext;
  Aext_mac(abs(Aext_mac) < 1e-8) = 0;
  Aext_mac_cell(:,2:size(Aext_mac,2)+1) = num2cell(Aext_mac);
  delete(strcat(directory_for_saving,"/supplementary_ext.xls"))
  writecell(Aext_mac_cell,strcat(directory_for_saving,"/supplementary_ext.xls"))

  % Kinetic parameters
  mets_not_included_in_Monod_kinetics = options_kinetic_identification.mets_not_included_in_Monod_kinetics;
  [mets_in_kinetics,locc] = ismember(mets_not_included_in_Monod_kinetics,data_cext.mets_meas_in_cext);
  Parameters_cell = cell(2+size(Amac,2),1+2*size(Kinetic_parameters.activation_parameters,2));
  Parameters_cell(1,1) = {"wmax"}; Parameters_cell(2,1) = {"wmax"};
  Parameters_cell(3:end,1) = num2cell(Kinetic_parameters.max_reaction_rate);
  for k = 1:size(Kinetic_parameters.activation_parameters,2)
    Parameters_cell(1,2*k) = data_cext.mets_meas_in_cext(k);
    Parameters_cell(1,2*k+1) = data_cext.mets_meas_in_cext(k);
    Parameters_cell(2,2*k) = {"Activation"};
    Parameters_cell(2,2*k+1) = {"Inhibition"};
    for j = 1:size(Amac,2)
      Parameters_cell(2+j,2*k) = num2cell(Kinetic_parameters.activation_parameters(j,k));
      Parameters_cell(2+j,2*k+1) = num2cell(Kinetic_parameters.inhibition_parameters(j,k));
    end
  end
  delete(strcat(directory_for_saving,"/parameters.xls"))
  writecell(Parameters_cell,strcat(directory_for_saving,"/parameters.xls"))
  fprintf("Done.\n") 
  pause(2)

  %% Saving various information about the model in a .mat file 
  fprintf("Saving the .mat file containing the model in the folder %s...\n",directory_for_saving)
  name_matlab_saving_file = strcat(directory_for_saving,"/",strcat(name_save_file,".mat"));
  delete(name_matlab_saving_file)
  save(name_matlab_saving_file,'q_all_models_train','q_all_models_predict','EFMs','data_stoich','Kinetic_parameters','options_kinetic_identification','options_EFM','options_data','Aext_mac','Amac','errors_CG_and_model','Amac','Parameters_cell');
  pause(2)

  %% Saving various information about the model in an excel file 
  fprintf("Saving the .xls file containing the model in the folder %s...\n",directory_for_saving)
  pause(2)
  name_excel_saving_file = strcat(directory_for_saving,"/",strcat(name_save_file,".xls"));
  index_in_sheet_3 = 2;
  delete(name_excel_saving_file)
  Parameters_cell(1,1) = {"Maximal rate constants"}; Parameters_cell(2,1) = {"Maximal rate constants"};
  writecell(Aext_mac_cell,name_excel_saving_file,'Sheet',1); % we save the macrosocpic stoichiometric matrix
  writecell(Parameters_cell,name_excel_saving_file,'Sheet',2); % we save the kinetic parameters
  writecell({" Macroscopic reactions "},name_excel_saving_file,'Sheet',3,'Range','A1:A1'); % we save the macroscopic reactions with all the extracellular metabolites and the reactions which constitute each EFM
  writecell({" Reactions composing each macroscopic reaction "},name_excel_saving_file,'Sheet',3,'Range','B1:B1');  
  for k = 1:size(Aext_mac,2)
    writecell(set_of_macro_reactions_ext(k),name_excel_saving_file,'Sheet',3,'Range',strcat("A",num2str(index_in_sheet_3),":","A",num2str(index_in_sheet_3)));  
    writecell(all_reactions_for_each_EFM(:,k),name_excel_saving_file,'Sheet',3,'Range',strcat("B",num2str(index_in_sheet_3),":","B",num2str(index_in_sheet_3-1+length(all_reactions_for_each_EFM(:,k)))));  
    index_in_sheet_3 = index_in_sheet_3 + length(all_reactions_for_each_EFM(:,k));
  end
  
  index_in_sheet_4 = 2;
  writecell({" Macroscopic reactions with only the measured metabolites "},name_excel_saving_file,'Sheet',4,'Range','A1:A1'); % we save the macroscopic reactions with only the MEASURED extracellular metabolites and the reactions which constitute each EFM  
  writecell({" Reactions composing each macroscopic reaction "},name_excel_saving_file,'Sheet',4,'Range','B1:B1');  
  for k = 1:size(Amac,2)
    writecell(set_of_macro_reactions_meas(k),name_excel_saving_file,'Sheet',4,'Range',strcat("A",num2str(index_in_sheet_4),":","A",num2str(index_in_sheet_4)));  
    writecell(all_reactions_for_each_EFM(:,k),name_excel_saving_file,'Sheet',4,'Range',strcat("B",num2str(index_in_sheet_4),":","B",num2str(index_in_sheet_4-1+length(all_reactions_for_each_EFM(:,k)))));  
    index_in_sheet_4 = index_in_sheet_4 + length(all_reactions_for_each_EFM(:,k));
  end
  
  writetable(Table_relative_error,name_excel_saving_file,'Sheet',5); % we save the modeling relative errors
  writetable(Table_absolute_error,name_excel_saving_file,'Sheet',6); % we save the modeling absolute errors
  
  % we save the cell specific rate computed from the kinetic model
  N = data_qext.N;
  nb_data_per_training_media = data_qext.nb_data_per_training_media;
  nb_data_per_prediction_media = data_qext.nb_data_per_prediction_media;
  ic = data_qext.ic;
  days_column = data_qext.days_column;
  cell_sheet7 = cell(3+size(Ameas,1),N);
  kl_index = 1;
  for k = 1:length(nb_data_per_training_media)
    nb_data_in_training_media_k = nb_data_per_training_media(k);
    ind_days_media_k = (ic == k);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_in_training_media_k
      cell_sheet7(1,kl_index) = {"Training"};
      cell_sheet7(2,kl_index) = media(index_training_media(k));
      cell_sheet7(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  for k = 1:length(nb_data_per_prediction_media)
    nb_data_per_prediction_media_k = nb_data_per_prediction_media(k);
    ind_days_media_k = (ic == k);
    all_days_for_media_k = days_column(ind_days_media_k);
    for j = 1:nb_data_per_prediction_media_k
      cell_sheet7(1,kl_index) = {"Prediction"};
      cell_sheet7(2,kl_index) = media(index_prediction_media(k));
      cell_sheet7(3,kl_index) = all_days_for_media_k(j);
      kl_index =  kl_index + 1;
    end
  end
  writecell({"Cell specific rates from kinetic model"},name_excel_saving_file,'Sheet',7,'Range','A1');
  writecell(cell_sheet7,name_excel_saving_file,'Sheet',7,'Range','B1');
  writecell(mets_meas_in_qext',name_excel_saving_file,'Sheet',7,'Range','A4');
  writecell(num2cell([q_all_models_train(:,:,2),q_all_models_predict(:,:,2)]),name_excel_saving_file,'Sheet',7,'Range','B4');

  clc
  fprintf("The identification is finished!\n")
  fprintf("All the saved files can be found in the folder %s\n",directory_for_saving)
end



function [macroscopic_reaction,all_reactions_for_each_EFM] = compute_macro_reactions(Amac,mets_in_Amac,A,mets,EFMs)

% This function computes the macrosocpic reaction involving the extracellular metabolites and each reaction involved in each computed EFM
% It identified the negative and positive stoichiometric coefficients in order to place them either on the left hand side (reacatnts) or the right hand side (products) of the reaction

% Computation of the macroscopic reactions
n_v = size(EFMs,1); % number of fluxes/reactions
macroscopic_reaction = cell(size(Amac,2),1); % Cell vector containing the macroscopic reactions
for j = 1:size(Amac,2)
    indj = ~(Amac(:,j) == 0); % we look for the non-zeros entries in the j-th column of Amac (i.e. the metabolites participating in the j-th reaction)
    stoichio_coeff = Amac(indj,j);
    names_mets_to_write_in_the_macro_reactions = mets_in_Amac(indj);
    all_negative_stoichio_coeff = stoichio_coeff(stoichio_coeff <0);
    names_mets_with_neg_stoichio_coeff = names_mets_to_write_in_the_macro_reactions(stoichio_coeff <0);
    all_positive_stoichio_coeff = stoichio_coeff(stoichio_coeff >0);
    names_pos = names_mets_to_write_in_the_macro_reactions(stoichio_coeff >0);

    left_hand_side_terms = [];
    for k = 1:length(all_negative_stoichio_coeff)
        if k ~= length(all_negative_stoichio_coeff)
            left_hand_side_terms =  strcat(left_hand_side_terms, strcat( num2str(abs(all_negative_stoichio_coeff(k))),{' '}, cell2mat(names_mets_with_neg_stoichio_coeff(k)),{' '},'+',{' '}));
        else
            left_hand_side_terms =  strcat(left_hand_side_terms, strcat( num2str(abs(all_negative_stoichio_coeff(k))),{' '}, cell2mat(names_mets_with_neg_stoichio_coeff(k)),{' '}));     
        end
    end

    right_hand_side_terms = [];
    for k = 1:length(all_positive_stoichio_coeff)

        if k ~= length(all_positive_stoichio_coeff)
            right_hand_side_terms =  strcat(right_hand_side_terms, strcat( num2str(all_positive_stoichio_coeff(k)),{' '}, cell2mat(names_pos(k)),{' '},'+',{' '}));
        else
            right_hand_side_terms =  strcat(right_hand_side_terms, strcat( num2str(all_positive_stoichio_coeff(k)),{' '}, cell2mat(names_pos(k))));     
        end
    end
    macroscopic_reaction(j) = strcat(left_hand_side_terms,'->', {' '}, right_hand_side_terms);
end

% Computation of the different reactions composing each EFM
% For this purpose, we look at all the reactions which are active in each EFM
% Then, as for the macroscopic rates, we write the reactions depending on the sign of the stoichiometric coefficients
EFMs(abs(EFMs) < 1e-8) = 0;
nb_active_reaction_per_EFMs = sum(abs(EFMs) > 0,1);
nb_max = max(nb_active_reaction_per_EFMs);
all_reactions_for_each_EFM = cell(nb_max,size(Amac,2));

for j = 1:size(Amac,2)
  nb_active_reactions_EFM_j = nb_active_reaction_per_EFMs(j);
  ind_ej_non_zeros = nonzeros(diag(abs(EFMs(:,j))> 0)*(1:n_v)');
  ej = EFMs(:,j);
  for kk = 1:nb_active_reactions_EFM_j
    ekj = sign(ej(ind_ej_non_zeros(kk)))*A(:,ind_ej_non_zeros(kk));
    indj = ~(ekj == 0);
    stoichio_coeff =ekj(indj);
    names_mets_to_write_in_the_macro_reactions = mets(indj);
    all_negative_stoichio_coeff = stoichio_coeff(stoichio_coeff <0);
    names_mets_with_neg_stoichio_coeff = names_mets_to_write_in_the_macro_reactions(stoichio_coeff <0);
    all_positive_stoichio_coeff = stoichio_coeff(stoichio_coeff >0);
    names_pos = names_mets_to_write_in_the_macro_reactions(stoichio_coeff >0);

    left_hand_side_terms = [];
    for k = 1:length(all_negative_stoichio_coeff)
        if k ~= length(all_negative_stoichio_coeff)
            left_hand_side_terms =  strcat(left_hand_side_terms, strcat( num2str(abs(all_negative_stoichio_coeff(k))),{' '}, cell2mat(names_mets_with_neg_stoichio_coeff(k)),{' '},'+',{' '}));
        else
            left_hand_side_terms =  strcat(left_hand_side_terms, strcat( num2str(abs(all_negative_stoichio_coeff(k))),{' '}, cell2mat(names_mets_with_neg_stoichio_coeff(k)),{' '}));     
        end
    end

    right_hand_side_terms = [];
    for k = 1:length(all_positive_stoichio_coeff)

        if k ~= length(all_positive_stoichio_coeff)
            right_hand_side_terms =  strcat(right_hand_side_terms, strcat( num2str(all_positive_stoichio_coeff(k)),{' '}, cell2mat(names_pos(k)),{' '},'+',{' '}));
        else
            right_hand_side_terms =  strcat(right_hand_side_terms, strcat( num2str(all_positive_stoichio_coeff(k)),{' '}, cell2mat(names_pos(k))));     
        end
    end
    all_reactions_for_each_EFM(kk,j) = strcat(num2str(ind_ej_non_zeros(kk)),{' : '},{' '},left_hand_side_terms,'->', {' '}, right_hand_side_terms);
  end
end

end
