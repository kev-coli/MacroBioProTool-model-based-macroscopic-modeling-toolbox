%%  Nominal optimization with a given metabolic network model
%   author : Mirko Pasquini
%
%   This function allows the user to solve the following nominal
%   optimization problem 
%
%   minimize{u,c}   -(F-q(end-1))*q(end)              (Harvest rate)          
%   subject to      g1: u - u_upper <= 0            (Medium bounds)
%                   g2: -u + u_lower <= 0
%                   g3: c - c_upper <= 0            (Concentration bounds)
%                   g4: c + c_upper <= 0            
%                   g5: q - q_upper <= 0            (Rate bounds)
%                   g6: -q + q_lower <= 0
%                   g7: a'u - b <= 0                (Solubility constraint)
%                   h : Xv*q - F*c + F[u;u_hat] = 0 (Mass Balance Eq.)
%
%   The function works as an interface to the function
%   solve_nominal_medium_optimization_explicit. 
%
%   optim_results = NOMINAL_OPTIMIZATION_GIVEN_METABOLIC_MODEL(model, initial_conditions, index_data, constraints)
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
%               index_data.involved_met_indices : vector of indexes,
%                   related to an extended stoichiometric matrix, for the
%                   metabolites of interest for the optimization (either to
%                   have constraints enforced or to optimize related
%                   quantities). This is not a mandatory field and only
%                   have use if an extended stoichiometric matrix is used
%                   (with unmeasured metabolites). If not specified, the
%                   involved_met_indices is set to
%                   length(initial_conditions.qext(:,1)), i.e. the number
%                   of rows in the stoichiometric matrix
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
%       max_iterations: specify the maximum number of iterations for the
%           optimization. If this is not specified the default number will
%           be 1e3.
%       objective_type: specify the type of objective to optimize. The
%           available options are:
%           'harvest' : maximize the predicted harvest of a certain 
%               metabolite of interest (whose index is rate_index)
%           'rate-min' : minimize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%           'rate-max' : maximize the predicted rate of a certain
%               metabolite of interest (whose index is rate_index)
%
%   @outputs
%       optim_results : structure containing information on the results
%           of the optimization.
%               optim_results.media_optimized : matrix in which the i-th
%                   column corresponds to the medium obtained by the
%                   optimization procedure when starting from the i-th 
%                   initial condition.
%               optim_results.cext_optimized : matrix in which the i-th
%                   column corresponds to the (predicted) vector of extracellular
%                   metabolite concentrations, obtained by the optimization
%                   procedure (and corresponding to the i-th optimized
%                   medium) when starting from the i-th initial condition.
%               optim_results.qext_optimized : matrix in which the i-th
%                   column corresponds to the (predicted) vector of uptake-secretion
%                   rates as a result of the optimization procedure when
%                   starting from the i-th initial condition
%               optim_results.exit_flag_optimized : vector in which the i-th 
%                   component is the exit flag of the optimization procedure
%                   when starting from the i-th initial condition. Please
%                   refer to fmincon documentation for a description of the
%                   exit flags, but in short if the exit flag is
%                   positive then the optimization can be considered
%                   succesful.
%               optim_results.first_order_optimality_optimized : vector in
%                   which the i-th component is the first-order optimality
%                   criterion obtained at the end of the optimization
%                   procedure when starting from the i-th initial
%                   condition. Please refer to fmincon documentation for a
%                   description of the first-order optimality criteria, but
%                   in short if this value is close to 0 we can consider
%                   that the optimization procedure converged to a feasible
%                   local optimizer. (any value below 1e-4 can be considered
%                   acceptable).
%               optim_results.objective_function_value_optimized : vector
%                   in which the i-th component correspond to the value of
%                   the objective function at the solution of the
%                   optimization procedure when starting form the i-th
%                   initial condition.
%               optim_results.lagrangian_mult : structure containing the 
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

function optim_results = nominal_optimization_given_metabolic_model(model, initial_conditions, index_data, constraints, display_settings, objective_type, max_iterations)
    
    if nargin < 7
        if nargin < 6
            objective_type = 'harvest';
        end
        max_iterations = 1e3;
    end

    if ~isfield(index_data, 'involved_met_indices')
        index_data.involved_met_indices = 1:length(initial_conditions.qext(:,1));
    end

    model.Amac = model.Amac(index_data.involved_met_indices,:);

    % 1. Load model

    model_reordered = reorder_model_for_optimization(model,index_data.decision_metabolite_indices,index_data.coupled_met_indices);

    % 2. Load starting conditions : media, concentrations, rates
    
    medium_temp = initial_conditions.medium;
    cext_temp = initial_conditions.cext;
    qext_temp = initial_conditions.qext;

    [medium_data, cext_data, qext_data] = reorder_data_for_optimization(medium_temp, cext_temp, qext_temp, index_data.decision_metabolite_indices, index_data.coupled_met_indices);

    % 3. Load medium bounds, concentration bounds, rate bounds and solubility
    % values

    medium_bounds = constraints.medium_bounds;
    nu = size(medium_bounds,1);
    rate_bounds = constraints.rate_bounds;
    concentration_bounds = constraints.concentration_bounds;
    % use the same reorder_data_for_optimization function to reorder the
    % indices in the bounds
    [~, concentration_bounds, rate_bounds] = reorder_data_for_optimization(medium_data, concentration_bounds, rate_bounds, index_data.decision_metabolite_indices, index_data.coupled_met_indices); % medium bounds are different and already refer only to decision variables

    a_solub = constraints.a_solub;
    b_solub = constraints.b_solub;

    % 4. Optimize for each starting condition and store results
    n_cond_opt = size(initial_conditions.medium,2);
    media_optimized = zeros(size(medium_bounds,1),n_cond_opt);
    cext_optimized = zeros(size(concentration_bounds,1),n_cond_opt);
    qext_optimized = zeros(size(rate_bounds,1),n_cond_opt);
    exit_flag_optimized = zeros(1,n_cond_opt);
    first_order_optimality_optimized = zeros(1,n_cond_opt);
    objective_function_value_optimized = zeros(1,n_cond_opt);
    lambda_optimized = cell(1,n_cond_opt);
    % with the re-ordering internal to this interface the position of rate
    % of interest and growth rate might be changed. This sets it right for
    % the optimization
    not_decision_metabolite_indices = setdiff(1:size(qext_optimized,1), index_data.decision_metabolite_indices);
    not_decision_metabolite_indices_coupled_met = setdiff(index_data.coupled_met_indices, index_data.decision_metabolite_indices);
    not_decision_metabolite_indices = colvec(not_decision_metabolite_indices);
    not_decision_metabolite_indices_coupled_met = colvec(not_decision_metabolite_indices_coupled_met);
    index_data.rate_index = find([index_data.decision_metabolite_indices; not_decision_metabolite_indices]==index_data.rate_index);
    index_data.growth_index = find([index_data.decision_metabolite_indices; not_decision_metabolite_indices]==index_data.growth_index);


    for i = 1:n_cond_opt
%       %== Uncomment to display optimization n. ==
%       disp(strcat('Optimizing condition n.',num2str(i)));
        % === Optimization with Harvest-Based Objective ===
        if nargin < 5
            [z,lambda,fval,exitflag,first_ord] = solve_nominal_medium_optimization_explicit(model_reordered,medium_data(index_data.decision_metabolite_indices,i),cext_data(:,i),...
                                                       medium_bounds,concentration_bounds,rate_bounds,a_solub,b_solub,index_data.coupled_met_indices,...
                                                       'none',medium_data(nu+1:end,i),index_data.rate_index,index_data.growth_index, max_iterations, objective_type);
        else
            [z,lambda,fval,exitflag,first_ord] = solve_nominal_medium_optimization_explicit(model_reordered,medium_data(index_data.decision_metabolite_indices,i),cext_data(:,i),...
                                                       medium_bounds,concentration_bounds,rate_bounds,a_solub,b_solub,index_data.coupled_met_indices,...
                                                       display_settings,medium_data(nu+1:end,i),index_data.rate_index,index_data.growth_index, max_iterations, objective_type);            
        end
        media_optimized(:,i) = z(1:nu);
        c_opt = z(nu+1:end);
        c0 = build_u_total(c_opt(1:nu),c_opt(nu+1:end),index_data); % NOTE: an improper use of build_u_total, but it does the trick right...
        [cext_optimized(:,i),qext_optimized(:,i)] = run_virtual_experiment(build_u_total(media_optimized(:,i),initial_conditions.medium(not_decision_metabolite_indices_coupled_met),index_data),model,c0,index_data.coupled_met_indices,'none');
        exit_flag_optimized(i) = exitflag;
        first_order_optimality_optimized(i) = first_ord;
        objective_function_value_optimized(i) = fval;
        lambda_optimized{i} = lambda;
    end

    % 5. Save results in adequate structure optim_results
    optim_results.media_optimized = media_optimized;
    optim_results.cext_optimized = cext_optimized;
    optim_results.qext_optimized = qext_optimized;
    optim_results.exit_flag_optimized = exit_flag_optimized;
    optim_results.first_order_optimality_optimized = first_order_optimality_optimized;
    optim_results.objective_function_value_optimized = objective_function_value_optimized;
    optim_results.lagrangian_mult = lambda_optimized;
end