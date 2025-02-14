%   [Kinetic_parameters,q_model_train,q_model_predict,w_model_train] = Step_3_Monod_kinetic_identification_of_macroscopic_rates(Macro_rates,data_cext,data_qext,data_stoich,EFMs,options_kinetic_identification,options_data,options_EFMs)  
%
%   This code performs the kinetic identification with three substeps (Steps 3.a, 3.b and 3.c).
%
%   @inputs:
%       - Macro_rates : matrix containing the computed macrosocpic rates from the column generation
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
%       - EFMs : EFM matrix
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
%
%   @outputs:
%       - Kinetic_parameters : structure containing the identified kinetic parameters
%               -> Kinetic_parameters.activation_parameters : matrix containing all the identified activation parameters. The rows correspond to the macro-rates and the column to the concentration (input)
%               -> Kinetic_parameters.inhibition_parameters : matrix containing all the identified inhibition parameters. The rows correspond to the macro-rates and the column to the concentration (input)
%               -> Kinetic_parameters.max_reaction_rate : vector containing all the identified maximal rate constants of each macroscopic rate.
%               -> Kinetic_parameters.lambda_l : optimal regularization parameter selected by the cross-validation. If no regularization was desired, then lambda_l is set to 0.
%       - q_model_train : cell specific rates computed from the kinetic model and evaluated with the training data
%       - q_model_predict : cell specific rates computed from the kinetic model and evaluated with the prediction data
%       - w_model_train : macroscopic rates computed from the kinetic model and evaluated with the training data

function [Kinetic_parameters,q_model_train,q_model_predict,w_model_train] = Step_3_Monod_kinetic_identification_of_macroscopic_rates(Macro_rates,data_cext,data_qext,data_stoich,EFMs,options_kinetic_identification)

  %% Initialization and loading of the necessary variables

  % Load the different options and data required for Step 3
  average_data_per_condition = options_kinetic_identification.average_data_per_condition; % contain the variable which specifies if the data of a same condition should be averaged.
  ic = data_cext.ic; % vector of indices which allow to identify the condition for each measurement
  ic_training = data_cext.ic_training; % vector of indices which allow to identify the TRAINING condition for each measurement
  mets_meas_in_cext = data_cext.mets_meas_in_cext; 
  mets_not_included_in_Monod_kinetics = options_kinetic_identification.mets_not_included_in_Monod_kinetics;
  regularization = options_kinetic_identification.regularization;
  lambda_grid = options_kinetic_identification.lambda_grid;
  number_multistart_points = options_kinetic_identification.number_multistart_points;
  number_max_of_iterations_EM = options_kinetic_identification.number_max_of_iterations_EM;
  burn_in_period_EM = options_kinetic_identification.burn_in_period_EM;
  number_samples_per_iteration_EM = options_kinetic_identification.number_samples_per_iteration_EM;
  number_trials_sampling_EM = options_kinetic_identification.number_trials_sampling_EM;
  perturbation_variance_EM = options_kinetic_identification.perturbation_variance_EM;
  index_training_media = data_qext.index_training_media;
  index_prediction_media = data_qext.index_prediction_media;
  Ameas = data_stoich.Ameas; 
  threshold_variation_neutral_effect = options_kinetic_identification.threshold_variation_neutral_effect;
  threshold_fit_activation_or_inhibition = options_kinetic_identification.threshold_fit_activation_or_inhibition;

  
  % Compute the macroscopic stoichiometric matrix
  Amac = Ameas*EFMs; 
  Amac(abs(Amac) < 1e-8) = 0; % remove negligible entries

  % Store the necessary options for the Bayesian estimation in the vector parameters_EM
  parameters_EM = zeros(5,1);
  parameters_EM(1,1) = number_max_of_iterations_EM;
  parameters_EM(2,1) = burn_in_period_EM;
  parameters_EM(3,1) = number_samples_per_iteration_EM;
  parameters_EM(4,1) = number_trials_sampling_EM;
  parameters_EM(5,1) = perturbation_variance_EM;

  % Loading the cell speciic rate, macro-rate and concentration data
  cext = data_cext.cext; % load the concentration data
  qext_train = data_qext.qext_train'; % load the cell specific rate TRAINING data
  [~,n_cext] = size(cext);
  [n_macro_rates,~] = size(Macro_rates); 
  cext_train = data_cext.cext_train; % load the concentration TRAINING data
  cext_predict = data_cext.cext_predict; % load the concentration PREDICTION data
  cext_train(cext_train < 1e-4) = 1e-4; % Replace zero concentration by a small value. This prevents numerical instability during modeling
  cext_predict(cext_predict < 1e-4) = 1e-4; % Replace zero concentration by a small value. This prevents numerical instability during modeling
  input = cext_train;
  input_predict = cext_predict;
  
  % Exclude the concentration data which are not included in the kinetics (decided by the user with the cell vector mets_not_included_in_Monod_kinetics,mets_meas_in_cext)
  [~,locc] = ismember(mets_not_included_in_Monod_kinetics,mets_meas_in_cext);
  locc = nonzeros(locc); ind_cext_in_kinetic_identification = setdiff(1:n_cext,locc);
  input = input(:,ind_cext_in_kinetic_identification); % input vector for the Monod model. Each column corresponds to a concentration of a metabolite which possibly participate in the macroscopic reactions
  n_inputs = size(input,2);
  output = Macro_rates;
  
  if(regularization == 1)
    %%%%%%%%%%%%%%%%%%%%%%  REGULARIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Leave-One-Out cross validation
    average_error_cross_valid = zeros(length(lambda_grid),1); % This vector contains the average testing error computed for each egularization parameter specified by the user in the variable lambda_grid 
    for l = 1:length(lambda_grid) % For each value of the regularization parameter specified by the user in the variable lambda_grid
      lambda_l = lambda_grid(l);
      for k = 1:length(index_training_media) % For each condition considered as test media (for the leave-one out cross validation)
        clc
        fprintf("Step 3 - Identification of the kinetics - Cross-validation with regularization.\n")
        fprintf("Regularization parameter %d/%d.\n",l,length(lambda_grid))
        fprintf("Validation medium %d/%d.\n",k,length(index_training_media))  
       
 	    % Split the training data between train and test data for the cross vdalidation
		ind_test_k = (ic_training == k); ind_train_k = logical(1 - ind_test_k); 
        input_test = input(ind_test_k,:); input_train = input(ind_train_k,:);
        output_train = Macro_rates(:,ind_train_k);
        qext_test_model_k = zeros(size(Amac,1),size(input_test,1));
        qext_test_k = qext_train(:,ind_test_k);
		
		% We identify a Monod model for each macroscopic reaction
		% Parallel computing can be used in order to accelerate the computation
        parfor j = 1:n_macro_rates 
          output_j_train = output_train(j,:)';
          if(sum(output_j_train<1e-7) < length(output_j_train))
            [~,~,~,~,w_test_model_j] = Monod_kinetic_identification(output_j_train,input_train,input_test,parameters_EM,number_multistart_points,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition,lambda_l);
            qext_test_model_k = qext_test_model_k + Amac(:,j)*w_test_model_j';
          end
        end
        average_error_cross_valid(l,1) = average_error_cross_valid(l,1) + norm(qext_test_model_k-qext_test_k,2)^2/size(qext_test_k,1);
      end
    end

    % Select regularization parameter with smallest average validation error
    % and run once again the identification. This will set many
    % parameters to zero. 
    activation_parameters = zeros(n_macro_rates,n_inputs);
    inhibition_parameters = zeros(n_macro_rates,n_inputs);
    max_reaction_rate = zeros(n_macro_rates,1);
    [~,ind_lambda_opt] = min(average_error_cross_valid);
    lambda_l = lambda_grid(ind_lambda_opt);
	
	% We identify a Monod model for each macroscopic reaction
    % Parallel computing can be used in order to accelerate the computation
    parfor j = 1:n_macro_rates
       output_j = output(j,:)'; 
       if(sum(output_j<1e-7) < length(output_j)) % if all the macroscopic rate data have negligible values, we do not fit a Model (the corresponding Monod parameters are all kept equal yo 0). Otherwize, we fit a Monod model
         [act_j,inh_j,alp_j,~,~] = Monod_kinetic_identification(output_j,input,[],parameters_EM,number_multistart_points,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition,lambda_l);
         activation_parameters(j,:) = act_j'; inhibition_parameters(j,:) = inh_j'; max_reaction_rate(j,1) = alp_j;
       end
    end

  else %%%%%%%%%%%%%%%%%%%% NO REGULARIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% In the case where the user does not want regularization, the code will simply try each condition as testing condition and compute the total error (training and test error)
    % The model with the smallest total error will be selected by the code
	lambda_l = 0;  
    activation_parameters = zeros(n_macro_rates,n_inputs);
    inhibition_parameters = zeros(n_macro_rates,n_inputs);
    max_reaction_rate = zeros(n_macro_rates,1);
    min_error_for_fitting_qext = Inf; 
    for k = 1:length(index_training_media)  % For each condition considered as test media 
      clc
      fprintf("Step 3 - Identification of the kinetics - Cross-validation without regularization.\n")
      fprintf("Validation medium %d/%d.\n",k,length(index_training_media))    
      % The data are split into trainign and testing data sets.
	  % The training set is used to compute the parameters_EM
	  % The test set is used to select the identified model
	  ind_test_k = (ic_training == k); ind_train_k = logical(1 - ind_test_k);
      input_test = input(ind_test_k,:); input_train = input(ind_train_k,:);
      output_test = output(:,ind_test_k); output_train = output(:,ind_train_k);
      activation_parameters_k = zeros(n_macro_rates,n_inputs);
      inhibition_parameters_k = zeros(n_macro_rates,n_inputs);
      max_reaction_rate_k = zeros(n_macro_rates,1);
      error_train_k = 0; error_test_k = 0;
      qext_train_model_k = zeros(size(Amac,1),size(input_train,1));
      qext_train_k = qext_train(:,ind_train_k); % Data used for model training
      qext_test_model_k = zeros(size(Amac,1),size(input_test,1));
      qext_test_k = qext_train(:,ind_test_k); % Data used for computing the testing error in order to guide model selection
	  
      % We identify a Monod model for each macroscopic reaction
	  % Parallel computing can be used in order to accelerate the computation
	  parfor j = 1:n_macro_rates
        output_j_test = output_test(j,:)'; output_j_train = output_train(j,:)';
        [act_jk,inh_jk,alp_jk,w_train_model_j,w_test_model_j] = Monod_kinetic_identification(output_j_train,input_train,input_test,parameters_EM,number_multistart_points,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition,lambda_l);
        activation_parameters_k(j,:) = act_jk'; inhibition_parameters_k(j,:) = inh_jk';  max_reaction_rate_k(j,1) = alp_jk;      
        qext_test_model_k = qext_test_model_k + Amac(:,j)*w_test_model_j';
        qext_train_model_k = qext_train_model_k + Amac(:,j)*w_train_model_j'; 
      end
	  % Computation of the total error (training error + test error)
      error_total_k = norm(qext_test_model_k-qext_test_k,2)^2/size(qext_test_k,1) + norm(qext_train_model_k-qext_train_k,2)^2/size(qext_train_k,1);
      if(error_total_k  < min_error_for_fitting_qext) % If the new model has a total error less than the previously best model, we store this model as best model
        min_error_for_fitting_qext = error_total_k;     
        best_valid_media = k;
        activation_parameters = activation_parameters_k;
        inhibition_parameters = inhibition_parameters_k;
        max_reaction_rate = max_reaction_rate_k;
      end 
    end

  end
  
  Kinetic_parameters.activation_parameters = activation_parameters;
  Kinetic_parameters.inhibition_parameters = inhibition_parameters;
  Kinetic_parameters.max_reaction_rate = max_reaction_rate;
  Kinetic_parameters.lambda_l = lambda_l;

  % Computation of the macroscopic rates computed from the kinetic model evaluated with the training and prediction data
  % w_model_train is the matrix of dimension containing the macroscopic rates computed using the kinetic model and evaluated with the training data
  % w_model_predict is the matrix containing the macroscopic rates computed using the kinetic model and evaluated with the prediction data
  % For both matrices, the rows correspond to the different macroscopic reactions and the columns to the different measurement time instants (data)
  N_train = size(input,1); N_predict = size(input_predict,1);
  w_model_train = ones(n_macro_rates,N_train); w_model_predict = ones(n_macro_rates,N_predict);
  for j = 1:n_macro_rates
     act_j = activation_parameters(j,:)'; inh_j = inhibition_parameters(j,:)'; alp_j = max_reaction_rate(j,1);
     for k = 1:n_inputs
       w_model_train(j,:) = w_model_train(j,:).*(input(:,k)')./(input(:,k)' + act_j(k))./(1 + inh_j(k)*input(:,k)');
       if(N_predict > 0)
         w_model_predict(j,:) = w_model_predict(j,:).*(input_predict(:,k)')./(input_predict(:,k)' + act_j(k))./(1 + inh_j(k)*input_predict(:,k)');
       end
     end
     w_model_train(j,:) = w_model_train(j,:)*alp_j;
     if(N_predict > 0)
       w_model_predict(j,:) = w_model_predict(j,:)*alp_j;
     end
  end
  
  % If average_data_per_condition is set to 1, the output of the model is the same for all the measurements of a same condition
  % The following part of the code just duplicates the value computed from the model for all the measurements of a same condition.
  % This is done for training and prediction conditions
  if(average_data_per_condition == 1)
    W_model_train = []; W_model_predict = [];
    kl_train = 1;
    for k = index_training_media
      ind_k = (ic == k); nb_data_in_condition_k = sum(ind_k);
      for l = 1:nb_data_in_condition_k
        W_model_train = [W_model_train,w_model_train(:,kl_train)];
      end
      kl_train = kl_train + 1;
    end
    kl_predict = 1;
    if(~isempty(index_prediction_media))
      for k = index_prediction_media
        ind_k = (ic == k); nb_data_in_condition_k = sum(ind_k);
        for l = 1:nb_data_in_condition_k
          W_model_predict = [W_model_predict,w_model_train(:,kl_predict)];
        end
        kl_predict = kl_predict + 1;
      end
    end
    w_model_train =  W_model_train; w_model_predict =  W_model_predict;
  end
   
  % q_model_train : Matrix containing the cell specific rates computed from the kinetic model and evaluated with the training data
  % the rows correspond to the metabolites and the columns to the different measurement time instants
  q_model_train = Amac*w_model_train;

  if(N_predict > 0) 
    q_model_predict = Amac*w_model_predict; % cell specific rates computed from the kinetic model and evaluated with the prediction data
  else
    q_model_predict = [];
  end

 
end
