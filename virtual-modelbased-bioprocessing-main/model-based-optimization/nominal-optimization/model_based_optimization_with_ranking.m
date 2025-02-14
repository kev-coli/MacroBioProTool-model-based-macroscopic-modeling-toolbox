%%  Solve nominal optimization problem from several initial conditions and
%%  associate a score to each solution based on objective function value, 
%%  distance from a given dataset and sensitivity to medium variability
%
%   author : Mirko Pasquini 
%
%   This function solves the nominal optimization problem starting from
%   different initial conditions and associates to each local solution a
%   score based on: objective function value, distance from dataset and
%   solution sensitivity. The user can decide the relative importance of
%   these three factors, through lambda_ranking. As the solution
%   sensitivity we will consider the standard deviation of the 
%   considered objective function value, with respect to medium variability.
%
%   [optimization_results, optimization_results_ranked] = 
%       MODEL_BASED_OPTIMIZATION_WITH_RANKING(model, ...
%                                             initial_conditions, ...
%                                             index_data, ...
%                                             constraints, ...
%                                             objective_type,...
%                                             distance_evaluation_vars, ...
%                                             sensitivity_evaluation_vars, ...
%                                             lambda_ranking, ...
%                                             verbose)                                                                    
%
%   @inputs:
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide" 
%           for further info.
%       initial_conditions : a structure containing the initial conditions
%           from which the optimization is started
%               initial_conditions.medium : matrix in which each column correspond
%                   to the medium of a different experimental condition;
%               initial_conditions.cext : matrix in which each column 
%                   correspond to the extracellular metabolite concentrations 
%                   of a different experimental condition;
%               initial_conditions.qext : matrix in which each column correspond 
%                   to the extracellular metabolite uptake/secretion rates of 
%                   a different experimental condition;
%       index_data : a structure containing all information related to 
%           important indexes for the optimization.
%               index_data.decision_metabolite_indices : vector of indexes 
%                   of the metabolites in the medium vector that has to be
%                   optimized (i.e. these are the main decision variables
%                   of the optimization problem). Please refer to the 
%                   "Beginner's guide" for a definition of decision metabolite.
%               index_data.coupled_met_indices : vector of indexes for the 
%                   coupled metabolites. The order refers to the order of
%                   metabolites in the Amac matrix. The usual case is that 
%                   the first (n-2) metabolites are coupled, leaving mAb and 
%                   biomass uncoupled.
%               index_data.rate_index : index of the rate of interest in the
%                   uptake/secretion rates vector
%               index_data.growth_index : index of the growth rate in the 
%                   uptake/secretion rates vector 
%       constraints : a structure containing all information related to the 
%           constraints of the optimization problem
%           constraints.medium_bounds : a matrix of dimension (nu x 2), 
%               where nu is the number of metabolites indexed by 
%               index_data.decision_metabolite_indices with the first column
%               being the lower bounds of the decision variables medium, 
%               and the second column being their upper bounds.
%           constraints.concentration_bounds : a matrix of dimension (nc x 2),
%               where nc is the number of extracellular metabolites indexed by
%               index_data.coupled_met_indices, with the first column being 
%               the lower bounds of those metabolites concentrations, and the
%               second column being their upper bounds.
%           constraints.rate_bounds : a matrix of dimension (nq x 2), where 
%               nq is the number of rates in the uptake/secretion rates
%               vector, with the first column being the lower bounds of 
%               those rates, and the second column being their upper bounds.
%           constraints.a_solub : linear coefficient of the solubility 
%               constraint a'*u <= b. Notice that u in this case refers only 
%               to the decision metabolites indexed by 
%               index_data.decision_metabolite_indices. 
%           constraints.b_solub  : constant term of the solubility constraint
%               a'*u <= b. Notice that u in this case refers only to the
%               decision metabolites indexed by index_data.decision_metabolite_indices.
%       objective_type : specify the type of objective to optimize. The
%           available options are:
%           'harvest' : maximize the predicted harvest of a certain 
%               metabolite of interest (whose index is rate_index)
%           'rate-min' : minimize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%           'rate-max' : maximize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%       distance_evaluation_vars : structure of variables used in the
%           distance evaluation of local solutions to a given dataset. 
%           Please refer to the function evaluate_distance_medium_from_dataset
%           for further information on distance evaluation of a medium from 
%           a given media dataset.
%               distance_evaluation_vars.norm_type : type of norm considered
%                   to evaluate the distance. The available options are the
%                   same as for the function 'norm' (see "help norm" for 
%                   further info), however common choices are '1' for the 
%                   1-norm and '2' for the 2-norm
%               distance_evaluation_vars.normalization_flag : if the flag is 
%                   set to 1, all media will be normalized so that their 
%                   components are between 0 and 1. This is done to reduce 
%                   the effect of having different interval lengths for 
%                   different metabolites
%               distance_evaluation_vars.media_data : matrix in which each 
%                   column corresponds to a medium in the data. Only the 
%                   components relative to decision variables of the 
%                   optimization problem are considered
%       sensitivity_evaluation_vars : structure of variables used in the
%           numerical sensitivity evaluation of local solutions. Please
%           refer to the function numerical_sensitivity_medium_variability 
%           for further information on distance evaluation of a medium from
%           a given media dataset.
%           sensitivity_evaluation_vars.met_names : metabolites labels 
%               from the stoichiometric matrix. The labels are in the form 
%               of a cell array of strings 
%               e.g. met_names = {'Glc'},{'Gln'},{'Ser'},...
%           sensitivity_evaluation_vars.uncertainty_range : a vector of 
%               length n_u, where n_u is the number of decision metabolites, 
%               containing the uncertainty range of the medium components with 
%               respect to which the variability of the concentrations, rates 
%               and harvest of the extracellular metabolites is evaluated. 
%               E.g.    uncertainty_range = [5,2,2,1], 
%                       nominal_medium = [55,15,10,8],
%               then the perturbed inputs will be generated in the range 
%                       [50, 60;
%                        13, 17;
%                        8, 12;
%                        7,  9];
%           sensitivity_evaluation_vars.N_realizations : number of perturbed
%               media generated in the range obtained by considering 
%               uncertainty_range. For each perturbed medium the 
%               corresponding steady-state concentrations and rates are 
%               evaluated, and based on their realizations, mean,
%               standard deviation, min and max are obtained.
%       lambda_ranking : structure containing the weights associated to
%           objective function value, distance from dataset and numerical
%           sensitivity, determining their importance in the evaluation of
%           the overall scores for the local solutions of the optimization
%               lambda_ranking.objective : weight associated to the
%                   objective function value of the optimization solutions
%               lambda_ranking.distance : weight associated to the distance
%                   of the optimization solutions from media_data 
%               lambda_ranking.sensitivity : weight associated to the
%                   sensitivity of the optimization solutions with respect
%                   to medium perturbation
%       verbose : 1 - verbose output active, 0 - verbose output inactive
%   @outputs:
%       optimization_results : structure containing all the information on
%           the local solutions obtained through the optimization
%           procedure, including information on the distance from the media
%           dataset and numerical sensitivity of each local solution. 
%               optimization_results.media_optimized : matrix in which 
%                   the i-th column corresponds to the medium obtained by 
%                   the optimization procedure when starting from the i-th 
%                   initial condition.
%               optimization_results.cext_optimized : matrix in which the 
%                   i-th column corresponds to the (predicted) vector of 
%                   extracellular metabolite concentrations, obtained by 
%                   the optimization procedure (and corresponding to the 
%                   i-th optimized medium) when starting from the i-th 
%                   initial condition.
%               optimization_results.qext_optimized : matrix in which the 
%                   i-th column corresponds to the (predicted) vector of 
%                   uptake-secretion rates as a result of the optimization 
%                   procedure when starting from the i-th initial condition
%               optimization_results.exit_flag_optimized : vector in which 
%                   the i-th component is the exit flag of the optimization 
%                   procedure when starting from the i-th initial condition.
%                   Please refer to fmincon documentation for a description 
%                   of the exit flags (but in short if the exit flag is
%                   positive then the optimization can be considered
%                   succesful).
%               optimization_results.first_order_optimality_optimized : 
%                   vector in which the i-th component is the first-order 
%                   optimality criterion obtained at the end of the 
%                   optimization procedure when starting from the i-th initial
%                   condition. Please refer to fmincon documentation for a
%                   description of the first-order optimality criteria, but
%                   in short if this value is close to 0 we can consider
%                   that the optimization procedure converged to a feasible
%                   local optimizer. (any value below 1e-4 can be considered
%                   acceptable).
%               optimization_results.objective_function_value_optimized : 
%                   vector
%                   in which the i-th component correspond to the value of
%                   the objective function at the solution of the
%                   optimization procedure when starting form the i-th
%                   initial condition.
%               optimization_results.lagrangian_mult : structure containing the 
%                   Lagrangian multipliers of the optimization
%                   lagrangian_mult.g1 : multipliers associated to the 
%                       upper bound constraints on the decision variable z,
%                       where z is the concatenation of u and c in the
%                       optimization problem
%                   lagrangian_mult.g2 : multipliers associated to the 
%                       lower bound constraints on the decision variable z
%                   lagrangian_mult.g3 : multipliers associated to the 
%                       upper bound constraints on the extracellular 
%                       metabolite rates q
%                   lagrangian_mult.g4 : multipliers associated to the 
%                       lower bound constraints on the extracellular 
%                       metabolite rates q
%                   lagrangian_mult.g5 : multipliers associated to the 
%                       solubility constraint
%                   lagrangian_mult.h : multipliers associated to the 
%                       mass-balance equation equality constraints
%               optimization_results.distance_scores : vector in which the
%                   i-th component is the distance of the i-th local
%                   solution (the i-th column of the matrix
%                   media_optimized) from the media dataset (expressed as
%                   the matrix media_data). Please refer to the function
%                   evaluate_distance_medium_from_dataset for further
%                   information.
%               optimization_results.index_nearest_media : vector in which
%                   i-th component is the index of the column in media_data
%                   which is closest to the i-th local solution of the
%                   optimization (i.e. the one that defines the distance)
%               optimization_results.sensitivity_scores : vector in which
%                   the i-th component is the score associated to the
%                   sensitivity of the i-th local solution of the 
%                   optimization. In particular the score is the
%                   sample-based standard deviation of the predicted
%                   objective function value corresponding to the i-th solution
%               optimization_results.sensitivity_structures : vector of
%                   cells in which the i-th entry is the sensitivity
%                   structure associated with the i-th local solution of
%                   the optimization. Please refer to the function
%                   numerical_sensitivity_medium_variability for further
%                   information of such structure.
%               optimization_results.solutions_score : vector in which the
%                   i-th component is the overall score associated to the
%                   i-th local solution of the optimization. This is
%                   evaluated as lambda_ranking.objective*optimization_scores + ...
%                   + lambda_ranking.distance*distance_scores ...
%                   - lambda_ranking.sensitivity*sensitivity_scores
%                   where optimization_scores is
%                   optimization_results.objective_function_value_optimized
%                   in the case of objective_type equal to 'harvest' or 
%                   'rate-max' and 
%                   -optimization_results.objective_function_value_optimized
%                   in the case of objective_type equal to 'rate-min'.
%       optimization_results_ranked : same as optimization_results but the
%           local solutions are ordered based on their score (in descending
%           order)


function [optimization_results, optimization_results_ranked] = model_based_optimization_with_ranking(model, ...
                                                                    initial_conditions, ...
                                                                    index_data, ...
                                                                    constraints, ...
                                                                    objective_type,...
                                                                    distance_evaluation_vars, ...
                                                                    sensitivity_evaluation_vars, ...
                                                                    lambda_ranking, ...
                                                                    verbose)

    if verbose == 1
        display_opt = 'iter-detailed';
    else
        display_opt = 'none';
    end

    % 1. Execute nominal optimization
    fprintf('Executing nominal optimization... \n');
    optimization_results = nominal_optimization_given_metabolic_model(model, initial_conditions, index_data, constraints, display_opt, objective_type);
    if strcmp(objective_type,'rate-min')
        optimization_scores = -optimization_results.objective_function_value_optimized; % in a case of rate minimization, to a higher rate correspond a lower score
    else % either rate or harvest maximization
        optimization_scores = optimization_results.objective_function_value_optimized;
    end
    % 2. For each solution evaluate distance from dataset
    fprintf('Evaluating distances from dataset... \n');
    [distance_scores, index_nearest_media] = evaluate_distance_scores(optimization_results, distance_evaluation_vars, constraints);

    % 3. For each solution evaluate sensitivity and return standard
    %    deviation of harvest
    fprintf('Evaluating sensitivity of local solutions... \n');
    [sensitivity_scores, sensitivity_structures] = evaluate_sensitivity_scores(optimization_results, sensitivity_evaluation_vars, model, index_data, initial_conditions, objective_type, verbose);

    % 4. Evaluate score for each solution
    solutions_score = lambda_ranking.objective*optimization_scores + ...
                      lambda_ranking.distance*distance_scores - ...
                      lambda_ranking.sensitivity*sensitivity_scores;

    [solutions_score_ranked, ranked_indexes] = sort(solutions_score, 'descend');

    optimization_results_ranked = rank_optimization_results(optimization_results, ranked_indexes);

    optimization_results_ranked.distance_scores = distance_scores(ranked_indexes);
    optimization_results.distance_scores = distance_scores;
    optimization_results_ranked.index_nearest_media = index_nearest_media(ranked_indexes);
    optimization_results.index_nearest_media = index_nearest_media;

    optimization_results_ranked.sensitivity_scores = sensitivity_scores(ranked_indexes);
    optimization_results.sensitivity_scores = sensitivity_scores;
    optimization_results_ranked.sensitivity_structures = sensitivity_structures(ranked_indexes);
    optimization_results.sensitivity_structures = sensitivity_structures;

    optimization_results_ranked.solutions_score = solutions_score_ranked; % they are already sorted
    optimization_results.solutions_score = solutions_score;
end

function [distance_scores, index_nearest_media] = evaluate_distance_scores(optimization_results, distance_evaluation_vars, constraints)
    n_conditions_optimized = size(optimization_results.media_optimized,2);
    distance_scores = zeros(1, n_conditions_optimized);
    index_nearest_media = distance_scores; % initialized as a vector with all zeros
    norm_type = distance_evaluation_vars.norm_type;
    normalization_flag = distance_evaluation_vars.normalization_flag;
    media_data = distance_evaluation_vars.media_data;
    medium_bounds = constraints.medium_bounds;
    
    for k = 1 : n_conditions_optimized
        medium = optimization_results.media_optimized(:,k);
        [distance_scores(k), index_nearest_media(k)] = evaluate_distance_medium_from_dataset(medium, media_data, norm_type, normalization_flag, medium_bounds);
    end
end

function [sensitivity_scores, sensitivity_structures] = evaluate_sensitivity_scores(optimization_results, sensitivity_evaluation_vars, model, index_data, initial_conditions, objective_type, verbose)
    n_conditions_optimized = size(optimization_results.media_optimized,2);
    met_names = sensitivity_evaluation_vars.met_names;
    uncertainty_range = sensitivity_evaluation_vars.uncertainty_range;
    N_realizations = sensitivity_evaluation_vars.N_realizations;
    sensitivity_scores = zeros(1, n_conditions_optimized);
    sensitivity_structures = cell(1, n_conditions_optimized);
    not_decision_metabolites = setdiff(index_data.coupled_met_indices, index_data.decision_metabolite_indices);
    u_hat = initial_conditions.medium(not_decision_metabolites,1); % extract fixed components of medium composition vector

    for k = 1 : n_conditions_optimized
        nominal_medium = optimization_results.media_optimized(:,k);
        nominal_concentration = optimization_results.cext_optimized(:,k);
        numerical_sensitivity_k = numerical_sensitivity_medium_variability(model, ...
                                                                met_names, ...
                                                                nominal_medium, ...
                                                                nominal_concentration, ...
                                                                uncertainty_range, ...
                                                                N_realizations, ...
                                                                u_hat, ...
                                                                index_data,...
                                                                verbose);
        sensitivity_structures{k} = numerical_sensitivity_k;
        if strcmp(objective_type,'harvest')
            sensitivity_scores(k) = numerical_sensitivity_k{index_data.rate_index}.harvest.std; % standard deviation of the harvest of the metabolite of interest
        else % either maximization or minimization of the uptake-secretion rate
            sensitivity_scores(k) = numerical_sensitivity_k{index_data.rate_index}.rates.std; % standard deviation of the uptake-secretion rate of the metabolite of interest
        end
    end
end

function optimization_results_ranked = rank_optimization_results(optimization_results, ranked_indexes)
    optimization_results_ranked.media_optimized = optimization_results.media_optimized(:,ranked_indexes);
    optimization_results_ranked.cext_optimized = optimization_results.cext_optimized(:,ranked_indexes);
    optimization_results_ranked.qext_optimized = optimization_results.qext_optimized(:,ranked_indexes);
    optimization_results_ranked.exit_flag_optimized = optimization_results.exit_flag_optimized(:,ranked_indexes);
    optimization_results_ranked.first_order_optimality_optimized = optimization_results.first_order_optimality_optimized(:,ranked_indexes);
    optimization_results_ranked.objective_function_value_optimized = optimization_results.objective_function_value_optimized(:,ranked_indexes);
    optimization_results_ranked.lagrangian_mult = optimization_results.lagrangian_mult(ranked_indexes);
end
