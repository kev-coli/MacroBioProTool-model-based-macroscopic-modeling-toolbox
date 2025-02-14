%% Rank optimization results based on lambda-ranking
%   author  :   Mirko Pasquini 
%
%   The function takes an optimization_results structure (for example as 
%   produced by the functions 
%   nominal_optimization_given_metabolic_model or 
%   model_based_optimization_with_ranking) and a lambda_ranking, describing
%   the relative importance of objective function value, distance from a given
%   dataset and sensitivity of the optimization solutions, and rank them
%   accordingly.
%   
%   This function can be used to re-rank the optimization results obtained
%   using the function model_based_optimization_with_ranking if we want to
%   consider a different lambda_ranking.
%	
%   optimization_results_ranked = GET_OPTIMIZATION_RESULTS_RANKING...
%                               (optimization_results, lambda_ranking)
%
%   @inputs:
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
%
%   @outputs:
%       optimization_results_ranked : same as optimization_results but the
%           local solutions are ordered based on the score obtained 
%           considering lambda_ranking (in descending order)


function optimization_results_ranked = get_optimization_results_ranking(optimization_results, lambda_ranking)
    
    optimization_scores = optimization_results.objective_function_value_optimized;
    distance_scores = optimization_results.distance_scores;
    sensitivity_scores = optimization_results.sensitivity_scores;

    solutions_score = lambda_ranking.objective*optimization_scores + ...
                      lambda_ranking.distance*distance_scores - ...
                      lambda_ranking.sensitivity*sensitivity_scores;

    [solutions_score_ranked, ranked_indexes] = sort(solutions_score, 'descend');

    optimization_results_ranked.media_optimized = optimization_results.media_optimized(:,ranked_indexes);
    optimization_results_ranked.cext_optimized = optimization_results.cext_optimized(:,ranked_indexes);
    optimization_results_ranked.qext_optimized = optimization_results.qext_optimized(:,ranked_indexes);
    optimization_results_ranked.exit_flag_optimized = optimization_results.exit_flag_optimized(:,ranked_indexes);
    optimization_results_ranked.first_order_optimality_optimized = optimization_results.first_order_optimality_optimized(:,ranked_indexes);
    optimization_results_ranked.objective_function_value_optimized = optimization_results.objective_function_value_optimized(:,ranked_indexes);
    optimization_results_ranked.lagrangian_mult = optimization_results.lagrangian_mult(ranked_indexes);
    optimization_results_ranked.distance_scores = optimization_results.distance_scores(ranked_indexes);
    optimization_results_ranked.index_nearest_media = optimization_results.index_nearest_media(ranked_indexes);
    optimization_results_ranked.sensitivity_scores = optimization_results.sensitivity_scores(ranked_indexes);
    optimization_results_ranked.sensitivity_structures = optimization_results.sensitivity_structures(ranked_indexes);
    optimization_results_ranked.solutions_score = solutions_score_ranked; % they are already sorted
end
