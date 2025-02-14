%% Numerical sensitivity with respect to medium variability
%   author  :   Mirko Pasquini 
%
%   The function is used to determine the sensitivity of some quantities
%   of interest (i.e. concentrations, rates and harvest of extracellular
%   metabolites) with respect to perturbations in a nominal medium 
%   composition. The sensitivity is evaluated numerically. In particular 
%   N_realizations perturbations of the nominal medium are generated, 
%   within the uncertainty_range, and min, max, mean and std are evaluated 
%   for all the quantities of interest.
%
%   numerical_sensitivity = NUMERICAL_SENSITIVITY_MEDIUM_VARIABILITY(model, ...
%                                                                   met_names, ...
%                                                                   nominal_medium, ...
%                                                                   nominal_concentration, ...
%                                                                   uncertainty_range, ...
%                                                                   N_realizations, ...
%                                                                   u_hat, ...
%                                                                   index_data,...
%                                                                   verbose)
% 
%   @input
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide"
%           for further info.
%       met_names : metabolites labels from the stoichiometric matrix. 
%           The labels are in the form of a cell array of strings 
%           e.g. met_names = {'Glc','Gln','Ser',,...
%       nominal_medium : vector of the extracellular metabolites 
%           concentrations in the medium, around which we are interested in
%           studying the concentrations, rates and harvest sensitivity. It
%           contains only the medium components relative to decision variables
%           of the optimization problem;
%       nominal_concentration : vector of extracellular metabolites
%           steady-state concentrations in the bioreactor associated to the
%           medium nominal_medium;
%       uncertainty_range : a vector of length n_u, where n_u =
%           length(nominal_medium), containing the uncertainty range of the 
%           medium components with respect to which the variability of the 
%           concentrations, rates and harvest of the extracellular metabolites
%           is evaluated. 
%           E.g.    uncertainty_range = [5,2,2,1], 
%                   nominal_medium = [55,15,10,8],
%           then the perturbed inputs will be generated in the range 
%                   [50, 60;
%                    13, 17;
%                     8, 12;
%                     7,  9];
%       N_realizations : number of perturbed media generated in the range
%           obtained by considering uncertainty_range. For each perturbed
%           medium the corresponding steady-state concentrations and rates
%           are evaluated, and based on all the realizations of these, mean,
%           standard deviation, min and max are evaluated.
%       u_hat : vector containing the components of the medium that are 
%           fixed, i.e. the ones that are not relative to decision variables 
%           of the optimization problem
%       index_data : a structure containing all information related to 
%           important indexes for the optimization.
%               index_data.decision_metabolite_indices : vector of indexes 
%                   of the metabolites in the medium vector that has to be
%                   optimized (i.e. these are the main decision variables
%                   of the optimization problem).
%               index_data.coupled_met_indices : vector of indexes for the 
%                   coupled metabolites. The order refers to the order of 
%                   metabolites in the Amac matrix. Please refer to the 
%                   "Beginner's guide" for a definition of coupled metabolite.
%               index_data.rate_index : index of the rate of interest in the
%                   uptake/secretion rates vector
%               index_data.growth_index : index of the growth rate in the 
%                   uptake/secretion rates vector 
%       verbose : 1 - output is printed with current simulation
%                 0 - no output is printed
%
%   @output
%       numerical_sensitivity : an array of structures containing information 
%           on how an uncertainty in the medium composition affects the 
%           concentrations, rates and harvest of the extracellular metabolites.
%           Each structure is associated to an extracellular metabolite, and 
%           it is composed of four fields:
%               name : contains the metabolite label
%               concentrations : contains the information relative to the
%                   variability of the steady-state concentration of the
%                   metabolite
%               rates : contains the information relative to the
%                   variability of the steady-state uptake-secretion rates 
%                   of the metabolite
%               harvest : contains the information relative to the
%                   variability of the steady-state harvest of the
%                   metabolite
%           Each of the above fields, with the exception of name, contains the
%           information relative to mean, standard deviation (std), min and
%           max of the considered quantity.
%           
%           REMARK : IT IS UP TO THE USER TO DETERMINE THE MEANING AND 
%               SIGNIFICANCE OF EACH QUANTITY. E.G. IT MIGHT NOT MAKE SENSE
%               TO CONSIDER THE HARVEST RATE OF GLUCOSE, SINCE THE HARVEST 
%               RATE IS A QUANTITY GENERALLY RELATED TO THE PRODUCT OF 
%               INTEREST YIELD.
%


function numerical_sensitivity = numerical_sensitivity_medium_variability(model, ...
                                                                met_names, ...
                                                                nominal_medium, ...
                                                                nominal_concentration, ...
                                                                uncertainty_range, ...
                                                                N_realizations, ...
                                                                u_hat, ...
                                                                index_data,...
                                                                verbose)
    n_c = length(nominal_concentration);
    nominal_rate = ComputeQExt(nominal_concentration, model.theta_matrix, model.Amac);
    n_q = length(nominal_rate);

    numerical_sensitivity = cell(max([n_c,n_q]),1);

    [~, concentration_realizations, rates_realizations] = ...
        run_montecarlo_simulations_around_medium(model, ...
                                                 nominal_medium, ...
                                                 nominal_concentration, ...
                                                 N_realizations,...
                                                 uncertainty_range, ...
                                                 u_hat, ...
                                                 index_data, ...
                                                 verbose);
    
    numerical_sensitivity = evaluate_sensitivities(numerical_sensitivity,...
                                                                concentration_realizations,...
                                                                rates_realizations,...
                                                                met_names,...
                                                                index_data,...
                                                                model);
                                                          
end

function [medium_realizations, concentration_realizations, rates_realizations] = ...
                        run_montecarlo_simulations_around_medium(   model,...
                                                                    nominal_medium, ...
                                                                    nominal_concentration, ...
                                                                    N_realizations,...
                                                                    uncertainty_range, ...
                                                                    u_hat, ...
                                                                    index_data, ...
                                                                    verbose)

    medium_realizations = zeros(length(nominal_medium),N_realizations);
    concentration_realizations = zeros(length(nominal_concentration), N_realizations);
    rates_realizations = zeros(size(model.Amac,1));

    for n = 1 : N_realizations
        if verbose == 1
            fprintf("Running simulation number %d of %d. \n",n,N_realizations);
        end
        medium_realizations(:,n) = nominal_medium - uncertainty_range + 2*uncertainty_range.*rand(size(uncertainty_range));
        total_medium = build_u_total(medium_realizations(:,n),u_hat,index_data);
        [concentration_realizations(:,n), rates_realizations(:,n)] = run_virtual_experiment(total_medium, model, nominal_concentration, index_data.coupled_met_indices, 'none');
    end
end

function numerical_sensitivity = evaluate_sensitivities(numerical_sensitivity,...
                                                        concentrations_realizations,...
                                                        rates_realizations,...
                                                        met_names,...
                                                        index_data,...
                                                        model)
    n_c = size(concentrations_realizations,1);
    n_q = size(rates_realizations,1);
    growth_rates_realizations = rates_realizations(index_data.growth_index,:);
    for k = 1 : max([n_c,n_q])
        if k <= n_c && k <= n_q
            numerical_sensitivity{k} = evaluate_quantity_sensitivity(concentrations_realizations(k,:), rates_realizations(k,:), growth_rates_realizations, met_names{k}, index_data, model);
        elseif k > n_c && k <= n_q
            numerical_sensitivity{k} = evaluate_quantity_sensitivity([], rates_realizations(k,:), growth_rates_realizations, met_names{k}, index_data, model);
        elseif k > n_q && k <= n_c % this will never happen I think 
            numerical_sensitivity{k} = evaluate_quantity_sensitivity(concentrations_realizations(k,:), growth_rates_realizations, [], met_names{k}, index_data, model);
        else
            fprintf('Error in the evaluation of the sensitivity. Index out of range. \n');
        end
    end
end

function sensitivity_k = evaluate_quantity_sensitivity(concentrations_realizations_k, rates_realizations_k, growth_rates_realizations, met_name_k, index_data, model)
    sensitivity_k.name = met_name_k;
    sensitivity_k.concentrations.min = min(concentrations_realizations_k);
    sensitivity_k.concentrations.max = max(concentrations_realizations_k);
    sensitivity_k.concentrations.mean = mean(concentrations_realizations_k);
    sensitivity_k.concentrations.std = std(concentrations_realizations_k);
    sensitivity_k.rates.min = min(rates_realizations_k);
    sensitivity_k.rates.max = max(rates_realizations_k);
    sensitivity_k.rates.mean = mean(rates_realizations_k);
    sensitivity_k.rates.std = std(rates_realizations_k);   
    harvest_rate_realizations = (model.F-growth_rates_realizations).*rates_realizations_k; % harvest rate of metabolite k
    sensitivity_k.harvest.min = min(harvest_rate_realizations);
    sensitivity_k.harvest.max = max(harvest_rate_realizations);
    sensitivity_k.harvest.mean = mean(harvest_rate_realizations);
    sensitivity_k.harvest.std = std(harvest_rate_realizations);
end
