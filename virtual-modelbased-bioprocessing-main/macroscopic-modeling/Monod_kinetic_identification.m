%   [act_final,inh_final,max_reaction_rate_final,w_train_model,w_test_model] = Monod_kinetic_identification(output_j_train,input_train,input_test,parameters_EM,number_multistart_points,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition,lambda_l)
%
%   This code contains the three subsetps Step 3.a, 3.b and 3.c. It is run for a specific macroscopic rate to be modeled as a Monod function.
%
%   - Step 3.a : it uses the Matlab nonlinear optimization function fmincon in order to minimized the regulqrized least-squares problem
%                It is combined with the Multistart function which allows to try different initial conditions for estimation improvement
%   - Step 3.b : it uses a custom made Matlab function called bayesian_estimation_Monod_kinetics.m
%   - Step 3.c : it uses a custom made Matlab function called Step_3_c_Monod_model_reduction.m
%
%   @inputs: 
%       - output_j_train : vector containing the data of the macroscopic rate to be modeled
%       - input_train : matrix containing the training concentration data. 
%       - input_test : matrix containing the testing concentration data (for model selection).
%       - parameters_EM : vector of 5 entries which specifies the different options chosen by the user for the EM algorithm
%       - number_multistart_points : number of initial parameter optimizations realized during Step 3.a with different initial conditions. The more initial identifications, the greater the chances for the computation of the global optimum
%       - threshold_variation_neutral_effect : threshold for reduction of an identified modulation function into a neutral effect. It should be non-negative and kept close to 0%
%                                              the closer to 0%, the less the chances to possibly reduce some identfiied functions into a neutral effect.
%       - threshold_fit_activation_or_inhibition : threshold for reduction of an identified modulation function into an activation or inhibition effect. It should be non-negative and kept close to 100% 
%                                                  the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%       - lambda_l : value of the regularization paremeter to be tried
%
%   @outputs: 
%       - act_final : vector of identified activation parameters. The j-th entry corresponds to the activation parameter related to the metabolite whose concentration is in the j-th column of input_train and input_test
%       - inh_final : vector of identified activation parameters. The j-th entry corresponds to the inhibition parameter related to the metabolite whose concentration is in the j-th column of input_train and input_test
%       - max_reaction_rate_final : identified maximal rate constant
%       - w_train_model : output of the Monod model evaluated with the training concentration
%       - w_test_model : output of the Monod model evaluated with the testing concentration (different from prediction data which are not involved in the model selection)

function [act_final,inh_final,max_reaction_rate_final,w_train_model,w_test_model] = Monod_kinetic_identification(output_j_train,input_train,input_test,parameters_EM,number_multistart_points,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition,lambda_l)
  
  warning('off')
  output_j_train_norm = output_j_train/max(output_j_train); % The macrosocpic rate is first normalized by dividing it by its maximal value
  n_input = size(input_train,2);
  
  % As initial point for the regularized least squares optimization, each parameter vector is chosen such that it is equal to the average of the corresponding concentration data
  act_init = mean(input_train',2); inh_init = act_init;
  w_init = ones(length(output_j_train_norm),1);
  for k = 1:n_input
    w_init = w_init.*input_train(:,k)./(input_train(:,k) + act_init(k))./(1 + input_train(:,k)*inh_init(k));
  end
  theta_init = [act_init;inh_init;w_init'*output_j_train_norm/(output_j_train_norm'*output_j_train_norm)];

  %% Step 3.a: model structure selection and initialization
  options_fmincon = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',10*10^3,'MaxFunctionEvaluations',10*10^3,'OptimalityTolerance',1e-8,'StepTolerance',1e-8,'FunctionTolerance',1e-8);
  ms = MultiStart; ms = MultiStart(ms,'Display','off','UseParallel',false);  
  problem = createOptimProblem('fmincon','objective',@(theta) least_squares_error_fit_Monod_kinetics(theta,output_j_train_norm,input_train,lambda_l),'x0',theta_init,'lb',[zeros(2*n_input,1);0],'ub',[10^3*max(input_train)';10^3*max(input_train)';Inf],'options',options_fmincon);  
  [theta_init_log,~] = run(ms,problem,number_multistart_points);
  act_init = theta_init_log(1:n_input); inh_init = theta_init_log((n_input+1):2*n_input);
  act_init(act_init < 1e-4*mean(input_train',2)) = 0; inh_init(inh_init < 1e-4*mean(input_train',2)) = 0; % Remove any negligible identified parameter

  %% Step 3.b: Bayesian estimation of the non-zero parameters of Step 3.a
  sample_act = (act_init > 0); % sample_act is a boleean vector used in Bayesian estimation. It indicates the activation parameters which should be tuned (1) and the ones which should be kept equal to 0 (0).
  sample_inh = (inh_init > 0); % sample_inh is a boleean vector used in Bayesian estimation. It indicates the inhibition parameters which should be tuned (1) and the ones which should be kept equal to 0 (0).
  [act_post,inh_post,~] = bayesian_estimation_Monod_kinetics(output_j_train_norm,input_train,parameters_EM,sample_act,sample_inh,act_init,inh_init);
  
  %% Step 3.c: Modulation function reduction
  [act_final,inh_final] = Step_3_c_Monod_model_reduction(input_train,act_post,inh_post,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition);
  
  %% Computation of the macroscopic rates from the final kinetic model
  N_train = size(input_train,1); N_test = size(input_test,1);
  w_train_model = ones(N_train,1); 
  for k = 1:n_input 
    w_train_model = w_train_model.*(input_train(:,k)./(input_train(:,k) + act_final(k))./(1 + inh_final(k)*input_train(:,k)));
  end
  max_reaction_rate_final = (w_train_model'*output_j_train)/(w_train_model'*w_train_model);
  w_train_model = max_reaction_rate_final*w_train_model;

  if(N_test >= 1)
    w_test_model = max_reaction_rate_final*ones(N_test,1);
    for k = 1:n_input 
      w_test_model = w_test_model.*(input_test(:,k)./(input_test(:,k) + act_final(k))./(1 + inh_final(k)*input_test(:,k)));
    end
  else
    w_test_model = [];   
  end

end


function [least_squares_error,gradient_least_squares_error] = least_squares_error_fit_Monod_kinetics(theta,output_train,input_train,lambda_l)

%   Function computing the regularized squared error between the macro rate data and the Monod model
%   It also computes the gradient of the regularized squared error with respect to the activation and inhibition parameters
%
%   @inputs: 
%       - theta : vector containing the activation and inhibition paremeters to be identified
%       - output_train : vector containing the training macroscopic rate data to be fitted with a Monod model
%       - input_train : matrix containing the training concentration data                                                the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%       - lambda_l : value of the regularization paremeter to be tried
%
%   @outputs: 
%       - least_squares_error : regularized square error between the Monod model computed with the parameter vector theta and the training and input training data
%       - gradient_least_squares_error : gradient of least_squares_error with respect to the parameters in theta

  n_input = size(input_train,2); % number of metabolites participating in the kinetics (number of inputs of the Monod model)
  act = theta(1:n_input); % vector containing all the activation parameters
  inh = theta(n_input+1:end-1); % vector containing all the inhibition parameters
  alp = theta(2*n_input+1,1);

  % Computation of the output of the Monod model with the training input data
  wmod0 = ones(length(output_train),1);
  for k = 1:n_input
    wmod0 = wmod0.*input_train(:,k)./(input_train(:,k) + act(k));
    wmod0 = wmod0./(1 + inh(k)*input_train(:,k));
  end
  wmod = alp*wmod0;
  
  % Computation of the squared error  between the output of the model and the training output data
  err = output_train - wmod;
  least_squares_error = norm(err,2)^2 + lambda_l*sum(act + inh);

  % Gradient of the least-squares cost with respect to the kinetic parameters
  gradient_least_squares_error = zeros(2*n_input+1,1);
  for k = 1:n_input
    gradient_least_squares_error(k,1) = 2*(wmod./(input_train(:,k) + act(k)))'*err + lambda_l;
    gradient_least_squares_error(n_input + k,1) = 2*(wmod.*input_train(:,k)./(1 + inh(k)*input_train(:,k)))'*err + lambda_l;
  end
  gradient_least_squares_error(2*n_input+1,1) = -2*wmod0'*err;

end


function [act_reduced,inh_reduced] = Step_3_c_Monod_model_reduction(input_train,act,inh,threshold_variation_neutral_effect,threshold_fit_activation_or_inhibition)
  
%   This function performs Step 3.c, i.e., modulation function reduction for a given Monod model (i.e., a given macroscopic rate)
%   For each modulation function of a given Monod model, it performs two tests
%   Test 1 : it first checks if the modulation function model can be reduced into a neutral effect
%   Test 2 : if Test 1 is unsucessful, a second test is performed in order to verify the possibility of a reduction into an activation or inhibition effect
%
%   @inputs: 
%       - input_train : matrix containing the training concentration data                                                the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%       - act : vector containing all the identified activation parameters of the corresponding Monod model
%       - inh : vector containing all the identified inhibition parameters of the corresponding Monod model
%       - threshold_variation_neutral_effect : threshold for reduction of an identified modulation function into a neutral effect. It should be non-negative and kept close to 0%
%                                              the closer to 0%, the less the chances to possibly reduce some identfiied functions into a neutral effect.
%       - threshold_fit_activation_or_inhibition : threshold for reduction of an identified modulation function into an activation or inhibition effect. It should be non-negative and kept close to 100% 
%                                                  the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%
%   @outputs: 
%       - act_reduced : vector containing all the identified activation parameters after modulation function reduction
%       - inh_reduced : vector containing all the identified inhibition parameters after modulation function reduction
  
  options = optimoptions('fmincon','SpecifyObjectiveGradient',false,'Display','none','MaxIterations',10^4,'MaxFunctionEvaluations',10^4,'OptimalityTolerance',1e-9,'StepTolerance',1e-9,'FunctionTolerance',1e-9);
  m = length(act);
  act_reduced = zeros(m,1); inh_reduced = zeros(m,1);
  for k = 1:m  % We perform the reduction for each modulation function

    % Because the double-component model is not identifiable, there are always two sets of parameters which fit the data either (act,inh) or (1/inh,1/act)
    % We select the one which has the smallest values for the activation and inhibition parameters	
    if(act(k)*inh(k) > 1)
      act_init_k = 1/inh(k);
      inh_init_k = 1/act(k);
    else
      act_init_k = act(k);
      inh_init_k = inh(k); 
    end

    xk = input_train(:,k); 
	
    % We evaluate the modulation function from the identified activation and inhibition parameters at the training inputs in xk	
    hk = xk./(xk + act_init_k)./(1 + inh_init_k*xk);

     % Test for neutral effect reduction
     % We check if the relative variation between the minimum and maximum in the vector hk is below an user-specified threshold
    if(100*(max(hk) - min(hk))/min(hk) < threshold_variation_neutral_effect) % If the test is successful, both activation and inhibition are set to 0
      act_reduced(k) = 0; inh_reduced(k) = 0; 
    else % If the test is not successful, we try to check if the modulation function can be reduced into an activation or inhibition effect
      
	  % We compute the activation function which fit the output of the modulation function model best (least squares optimization)
      act_k = fmincon(@(act) cost_least_squares_fit_modulation_function(act,xk,hk),act_init_k,[],[],[],[],0,[],[],options);
      hact0_k = xk./(xk + act_init_k); alp_act = (hact0_k'*hk)/(hact0_k'*hact0_k); 
      fit_act = 100*(1-norm(hk - alp_act*hact0_k,2)/norm(hk-mean(hk),2)); % Data fit of the activation function

      % We compute the inhibition function which fit the output of the modulation function model best (least squares optimization)
      inh_k = fmincon(@(inh) cost_least_squares_fit_modulation_function(inh,1./xk,hk),inh_init_k,[],[],[],[],0,[],[],options);
      hinh0_k = 1./(1 + inh_init_k*xk); alp_inh = (hinh0_k'*hk)/(hinh0_k'*hinh0_k); 
      fit_inh = 100*(1-norm(hk - alp_inh*hinh0_k,2)/norm(hk-mean(hk),2)); % Data fit of the inhibition function
      
      [fit_max,ind_max] = max([fit_act;fit_inh]); % We select the maximal fit between both models
	  
      if(fit_max < threshold_fit_activation_or_inhibition) % If the maximal fit is NOT above an user-specified threshold, then the modulation function CANNOT be reduced into the corresponding effect (activation or inhibition)
        act_reduced(k) = act_init_k; inh_reduced(k) = inh_init_k; 
      else % If the maximal fit is above an user-specified threshold, then the modulation function can be reduced into the corresponding effect (activation or inhibition)
        if(ind_max == 1) % If the maximal fit comes from the activation model
          act_reduced(k) = act_k; inh_reduced(k) = 0; 
        else % If the maximal fit comes from the inhibition model
          act_reduced(k) = 0; inh_reduced(k) = inh_k; 
        end
      end
    end
  end

end


function [cost] = cost_least_squares_fit_modulation_function(theta,x,h)

% theta: activation or inhibition parameter
% x: input, either c for activation or 1/c for inhibition (c = concentration)
% h: the modulation function data
%
%   @inputs: 
%       - theta : vector containing the activation and inhibition paremeters to be identified
%       - output_train : vector containing the training macroscopic rate data to be fitted with a Monod model
%       - input_train : matrix containing the training concentration data                                                the closer to 100%, the less the chances to possibly reduce some identified functions into an activation or inhibition effect.
%       - lambda_l : value of the regularization paremeter to be tried
%
%   @outputs: 
%       - least_squares_error : regularized square error between the Monod model computed with the parameter vector theta and the training and input training data
%       - gradient_least_squares_error : gradient of least_squares_error with respect to the parameters in theta
 
  hmod0 = x./(x+theta); 
  alp = h'*hmod0/(hmod0'*hmod0);
  cost = norm(h - alp*hmod0,2)^2;

end