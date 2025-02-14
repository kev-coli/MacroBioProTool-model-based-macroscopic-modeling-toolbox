%   [act_bayes,inh_bayes,max_reaction_rate] = bayesian_estimation_Monod_kinetics(output_j_train,input_train,parameters_EM,sample_act,sample_inh,act_init,inh_init) 
%
%   This code contains the Bayesian estimation algorithm (Step 3.b)
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

function [act_bayes,inh_bayes,max_reaction_rate] = bayesian_estimation_Monod_kinetics(output_j_train,input_train,parameters_EM,sample_act,sample_inh,act_init,inh_init) 
  
    % Options of the Expectation Maximization algorithm specified by the user
    number_max_of_iterations_EM = parameters_EM(1,1);
    burn_in_period_EM = parameters_EM(2,1);
    number_samples_per_iteration_EM = parameters_EM(3,1);
    number_trials_sampling_EM = parameters_EM(4,1);
    perturbation_variance_EM = parameters_EM(5,1);
      
    N = length(output_j_train); % number of data
    n_input = size(input_train,2); % number of metabolites involved in the macroscopic rate (inputs of the Monod model)
 
    % Initial values of all the hyperparameters (mean and standard deviations of the log-normal distribution for each activation and inhibition parameter)
    % The means are initialized by computing the log10 of the initial value of the activation and inhibition parameters
    act_log_mean = log10(act_init); % act_log_mean : vector containing the means of the log prior distribution of each activation parameter
    inh_log_mean = log10(inh_init); % inh_log_mean : vector containing the means of the log prior distribution of each inhibition parameter
    act_log_standard_dev = 1/2*ones(n_input,1); % act_log_standard_dev : vector containing the standard deviations of the log prior distribution of each activation parameter
    inh_log_standard_dev = 1/2*ones(n_input,1); % inh_log_standard_dev : vector containing the standard deviations of the log prior distribution of each inhibition parameter

    act_max = 1000*max(input_train)'; % Maximal value allowed for the sampling for the activation parameters
    inh_max = 1000*max(input_train)'; % Maximal value allowed for the sampling for the inhibition parameters
    
    theta_current_sample = zeros(2*n_input,1); % Parameter vector which gathers the recent sample of each activation an inhibition parameter
    % This parameter vector is initialized with the initial parameter vector coming from Step 3.a
    for kk = 1:n_input
      if(sample_act(kk,1) == 1)  
        theta_current_sample(kk,1) = act_init(kk);
      end
      if(sample_inh(kk,1) == 1)  
        theta_current_sample(n_input + kk,1) = inh_init(kk);
      end
    end
    theta_posterior = theta_current_sample; % theta_posterior : posterior estimate of the kinetic parameters. 
    act_posterior = theta_posterior(1:n_input);
    inh_posterior = theta_posterior(n_input+1:2*n_input);
    monod_model = ones(N,1);
    modulation_function = zeros(N,2*n_input);
    for kk = 1:n_input
      input_k = input_train(:,kk);  
      monod_model = monod_model.*input_k./(input_k+act_posterior(kk))./(1+inh_posterior(kk)*input_k);  
      modulation_function(:,kk) = input_k./(input_k+act_posterior(kk));
      modulation_function(:,n_input + kk) = 1./(1+inh_posterior(kk)*input_k);
    end
                
    max_reaction_rate = (monod_model'*output_j_train)/(monod_model'*monod_model);
    current_fit = 100*(1-norm(max_reaction_rate*monod_model - output_j_train,2)/norm(output_j_train - mean(output_j_train),2)); % Data fit of the new posterior estimate of the kinetic parameters
    best_fit = current_fit;
    noise_variance = 1/N*norm(max_reaction_rate*monod_model- output_j_train,2)^2;
    
    % At each iteration, the algorithm computes the data fit of the Monod model computed with the last posterior estimate of the kinetics parameters 
    % It stores the posterior estimate with the highest data fit in the variable best_theta_posterior and best_max_reaction_rate
    best_theta_posterior = theta_posterior;
    best_max_reaction_rate = max_reaction_rate;
    consecutive_failure_for_improving_fit = 0;

    % Kernel construction
    for jj = 1:number_max_of_iterations_EM
        
        if(jj > 1) % The burn-in period is done during the first iteration of the EM algorithm. It is set to 0 for the next iterations.
          burn_in_period_EM = 0;
        end

        %% E-step: sample kinetic parameters
        all_parameter_samples = zeros(2*n_input,number_samples_per_iteration_EM);
      
        for mm = (-burn_in_period_EM+1):number_samples_per_iteration_EM
          for kk = 1:n_input
            if(sample_act(kk) == 1) % We sample act first
              parameter_log_mean = act_log_mean(kk);
              parameter_log_standard_dev = act_log_standard_dev(kk);    
              parameter_previous_sample = theta_current_sample(kk,1);
              condition_set = setdiff(1:2*n_input,kk);
              product_other_terms = prod(modulation_function(:,condition_set),2);
              input_k = input_train(:,kk);
              theta_max = act_max(kk);
              [new_modulation_function,parameter_new_sample] = sample_parameter_metropolis_hastings(input_k, output_j_train, product_other_terms, 1, parameter_log_mean, parameter_log_standard_dev, parameter_previous_sample, noise_variance, number_trials_sampling_EM, theta_max); 
              modulation_function(:,kk) = new_modulation_function/mean(new_modulation_function);
              theta_current_sample(kk,1) = parameter_new_sample;
            end

            if(sample_inh(kk) == 1)  % We sample inh afterwards
              parameter_log_mean = inh_log_mean(kk);
              parameter_log_standard_dev = inh_log_standard_dev(kk);    
              parameter_previous_sample = theta_current_sample(n_input + kk,1);
              condition_set = setdiff(1:2*n_input,n_input + kk);
              product_other_terms = prod(modulation_function(:,condition_set),2);
              input_k = input_train(:,kk);
              theta_max = inh_max(kk);
              [new_modulation_function,parameter_new_sample] = sample_parameter_metropolis_hastings(input_k, output_j_train, product_other_terms, 2, parameter_log_mean, parameter_log_standard_dev, parameter_previous_sample, noise_variance, number_trials_sampling_EM, theta_max); 
              modulation_function(:,n_input + kk) = new_modulation_function/mean(new_modulation_function);
              theta_current_sample(n_input + kk,1) = parameter_new_sample;
            end     
          end
          if mm >0
            all_parameter_samples(:,mm) = theta_current_sample;
          end
        end

        %% Compute posterior model
        theta_posterior = mean(all_parameter_samples,2); 
        act_mean_post = theta_posterior(1:n_input,1);
        inh_mean_post = theta_posterior((n_input + 1):2*n_input,1);
        monod_model_post = ones(N,1);
        modulation_function_post = zeros(N,2*n_input);
        for kk = 1:n_input
          input_k = input_train(:,kk);  
          monod_model_post = monod_model_post.*input_k./(input_k+act_mean_post(kk))./(1+inh_mean_post(kk)*input_k);  
          modulation_function_post(:,kk) = input_k./(input_k+act_mean_post(kk));
          modulation_function_post(:,n_input + kk) = 1./(1+inh_mean_post(kk)*input_k);
        end
        
        % adapt the maximal constant rate
        max_reaction_rate = (monod_model_post'*output_j_train)/(monod_model_post'*monod_model_post);
        
        % new fit
        current_fit = 100*(1-norm(max_reaction_rate*monod_model_post - output_j_train,2)/norm(output_j_train - mean(output_j_train),2));
        
        % new estimate of the noise variance
        noise_variance = norm(output_j_train-max_reaction_rate*prod(monod_model_post,2),2)^2/N;
        
        if(current_fit > best_fit) % If the new posterior model has a better fit than the previously best model; then we store the information fo this new model
          best_fit = current_fit;
          best_max_reaction_rate = max_reaction_rate;
          best_theta_posterior = theta_posterior;
          consecutive_failure_for_improving_fit = 0;
        else
          consecutive_failure_for_improving_fit = consecutive_failure_for_improving_fit + 1;  
        end
     
        %% M-step: update hyperparameters from the samples of the kinetic parameters
        for kk=1:n_input
          if(sample_act(kk,1) == 1)
            act_log_mean(kk) = mean(log(all_parameter_samples(kk,:))/log(10)); % compute the mean of the log-distribution of each activation parameter
            act_log_standard_dev(kk) = sqrt(var(log(all_parameter_samples(kk,:))/log(10))) + perturbation_variance_EM; % compute the standard deviation of the log-distribution of each activation parameter
          end                                                                                                          % an additive term perturbation_variance_EM is added in order to explore the parameter space
          if(sample_inh(kk,1) == 1) 
            inh_log_mean(kk) = mean(log(all_parameter_samples(n_input + kk,:))/log(10)); % compute the mean of the log-distribution of each inhibition parameter
            inh_log_standard_dev(kk) = sqrt(var(log(all_parameter_samples(n_input + kk,:))/log(10))) + perturbation_variance_EM; % compute the standard deviation of the log-distribution of each inhibition parameter
          end                                                                                                                    % an additive term perturbation_variance_EM is added in order to explore the parameter space  
        end       
    end
    % At the end of the algorithm, the posterior estimate which is kept is the one with the maximal data fit
    max_reaction_rate = best_max_reaction_rate;
    theta_posterior = best_theta_posterior;  
    act_bayes = theta_posterior(1:n_input,1);
    inh_bayes = theta_posterior((n_input+1):2*n_input,1);
end



function [new_modulation_function,parameter_new_sample] = sample_parameter_metropolis_hastings(input_k, output, product_other_terms, type_of_kinetics, parameter_log_mean, parameter_log_standard_dev, parameter_previous_sample, noise_variance, number_trials_sampling_EM, theta_max)

%  This code handles the sampling used during the Bayesian estimation of the Monod parameters
%  It is based on the Metropolis-Hastings scheme  
%
%   @inputs: 
%       - input_k : vector containing the concentration training data of the corresponding metabolite whose kinetic parameter is sampled when this function is called
%       - output : vector of the macroscopic rate training data 
%       - product_other_terms : vector corresponding to the product of all the modulation functions except the one whose parameter is sampled
%       - type_of_kinetics : scalar defining the type of parameter sampled. 1 = activation, 2 = inhibition
%       - parameter_log_mean : 
%       - parameter_log_standard_dev : 
%       - parameter_previous_sample :
%       - noise_variance :
%       - number_trials_sampling_EM :
%       - theta_max :
%   @outputs: 
%       - new_modulation_function : vector of identified activation parameters. The j-th entry corresponds to the activation parameter related to the metabolite whose concentration is in the j-th column of input_train and input_test
%       - parameter_new_sample : vector of identified activation parameters. The j-th entry corresponds to the inhibition parameter related to the metabolite whose concentration is in the j-th column of input_train and input_test

previous_sample = parameter_previous_sample;
accept_sample = 0; % boolan variable. If it is equal to 1, the candidate sample has been accepted. Otherwise, the algorithm has not been able to sample one
n_attempts_at_sampling = 0; % number of current attempts at sampling an adequate parameter
log_previous_sample = log10(previous_sample);  % log10 of the previous sample of the kinetic parameter to be sampled
log_prior_previous_sample = -1/2*norm(log_previous_sample - parameter_log_mean,2)^2*noise_variance; % log prior evaluated at the previous sample
  
if(type_of_kinetics == 1) % if the parameter to be sampled is an activation parameter
  previous_modulation_function = diag(product_other_terms)*(input_k./(input_k+previous_sample)); %  modulation function evaluated with the previous sample
  alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);  
  log_likelihood_previous_sample = -1/2*norm(output-alp_h*previous_modulation_function,2)^2*parameter_log_standard_dev^2; % log likelihood evaluated at the previous sample
elseif(type_of_kinetics == 2) % if the parameter to be sampled is an inhibition parameter
  previous_modulation_function = diag(product_other_terms)*(1./(1+previous_sample*input_k)); %  modulation function evaluated with the previous sample
  alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);  
  log_likelihood_previous_sample = -1/2*norm(output-alp_h*previous_modulation_function,2)^2*parameter_log_standard_dev^2;  % log likelihood evaluated at the previous sample 
end
     
while(accept_sample == 0 && n_attempts_at_sampling < number_trials_sampling_EM)
  log_candidate_sample = log_previous_sample + parameter_log_standard_dev*randn(1,1); % log10 of the candidate sample for the parameter
  parameter_candidate = 10^(log_candidate_sample); % value of candidate sample for the parameter   
  if(type_of_kinetics == 1) % if the parameter is an activation parameter
    previous_modulation_function = diag(product_other_terms)*(input_k./(input_k+parameter_candidate));
    alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);
    log_likelihood_candidate_sample = - 1/2*norm(output-alp_h*previous_modulation_function,2)^2*parameter_log_standard_dev^2;
    log_prior_candidate_sample = - 1/2*norm(log_candidate_sample - parameter_log_mean,2)^2*noise_variance;    
    ratio_for_metropolis_hastings_test = log_likelihood_candidate_sample + log_prior_candidate_sample - log_likelihood_previous_sample - log_prior_previous_sample;
    sample_for_metropolis_hastings_test = rand(1,1); % sample for the Metropolis-Hastings test for candidate sample acceptance
    if(ratio_for_metropolis_hastings_test > log(sample_for_metropolis_hastings_test)*noise_variance*parameter_log_standard_dev^2)  % if the test is passed
      accept_sample = 1;
      new_modulation_function = alp_h*input_k./(input_k+parameter_candidate);
      previous_sample = parameter_candidate;
    else % if it fails
      n_attempts_at_sampling = n_attempts_at_sampling + 1;  
    end
  elseif(type_of_kinetics == 2) % if the parameter is an inhibition parameter
    previous_modulation_function = diag(product_other_terms)*(1./(1+parameter_candidate*input_k));
    alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);
    log_likelihood_candidate_sample = - 1/2*norm(output-alp_h*previous_modulation_function,2)^2*parameter_log_standard_dev^2;
    log_prior_candidate_sample = - 1/2*norm(log_candidate_sample - parameter_log_mean,2)^2*noise_variance;
    ratio_for_metropolis_hastings_test = log_likelihood_candidate_sample + log_prior_candidate_sample - log_likelihood_previous_sample - log_prior_previous_sample;
    sample_for_metropolis_hastings_test = rand(1,1); % sample for the Metropolis-Hastings test for candidate sample acceptance
    if(ratio_for_metropolis_hastings_test > log(sample_for_metropolis_hastings_test)*noise_variance*parameter_log_standard_dev^2)
      accept_sample = 1;
      new_modulation_function = alp_h*(1./(1+parameter_candidate*input_k));
      previous_sample = parameter_candidate;
    else
      n_attempts_at_sampling = n_attempts_at_sampling + 1;  
    end
  end 
end

parameter_new_sample = previous_sample;
if(parameter_new_sample > theta_max) % It happens that large parameters can be sampled. It can cause numerical issues if taken too large.
  parameter_new_sample = theta_max;
end

if(type_of_kinetics == 1)
  previous_modulation_function = diag(product_other_terms)*(input_k./(input_k+parameter_new_sample));
  alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);  
  new_modulation_function = alp_h*input_k./(input_k+parameter_new_sample);
elseif(type_of_kinetics == 2)
  previous_modulation_function = diag(product_other_terms)*(1./(1+parameter_new_sample*input_k));
  alp_h = (previous_modulation_function'*output)/(previous_modulation_function'*previous_modulation_function);  
  new_modulation_function = alp_h*1./(1+parameter_new_sample*input_k);  
end

end




