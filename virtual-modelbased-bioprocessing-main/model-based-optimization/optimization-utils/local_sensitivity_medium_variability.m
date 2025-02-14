%% Local sensitivity with respect to medium variability
%   author  :   Mirko Pasquini 
%
%   The function is used to determine the sensitivity of some quantities
%   of interest (i.e. concentrations, rates and harvest of extracellular
%   metabolites) with respect to perturbations in a nominal medium 
%   composition. The sensitivity is evaluated locally through the use of 
%   the Jacobian of the above quantities with respect to the medium 
%   components, evaluated at the nominal medium.
%
%   local_sensitivity = LOCAL_SENSITIVITY_MEDIUM_VARIABILITY(  model, ...
%                                                              met_names, ...
%                                                              nominal_medium, ...
%                                                              nominal_concentration, ...
%                                                              index_data)
%
%   @input
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide"
%           for further info.
%       met_names : metabolites labels from the stoichiometric matrix. 
%           The labels are in the form of a cell array of strings 
%           e.g. met_names = {'Glc','Gln','Ser',...
%       nominal_medium : vector of the extracellular metabolites 
%           concentrations in the medium, around which we are interested in
%           studying the sensitivity of the quantities of interest. It
%           contains only the medium components relative to decision variables
%           of the optimization problem;
%       nominal_concentration : vector of extracellular metabolites
%           steady-state concentrations in the bioreactor associated to the
%           medium nominal_medium;
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
%
%   @output
%       local_sensitivity : an array of structures containing information 
%           on how a variation in the medium composition affects the 
%           concentrations, rates and harvest of the extracellular metabolites.
%           Each structure is associated to an extracellular metabolite, and 
%           it iscomposed of four fields:
%               name : contains the metabolite label
%               concentrations : contains the information relative to the
%                   sensitivity of the steady-state concentration of the
%                   metabolite, with respect to a small variation in the
%                   medium composition;
%               rates : contains the information relative to the
%                   sensitivity of the steady-state uptake-secretion rates 
%                   of the metabolite, with respect to a small variation in
%                   the medium composition;
%               harvest : contains the information relative to the
%                   sensitivity of the steady-state harvest of the
%                   metabolite, with respect to a small variation in the
%                   medium composition.
%           Each of the above fields, with the exception of name, contains the
%           information relative local sensitivity of the considered
%           quantity with respect to the medium composition (i.e. it is the
%           gradient of the considered quantity with respect to the
%           components of medium relative to decision variables of the
%           optimization problem).
%           
%           REMARK : IT IS UP TO THE USER TO DETERMINE THE MEANING AND 
%               SIGNIFICANCE OF EACH QUANTITY. E.G. IT MIGHT NOT MAKE SENSE 
%               TO CONSIDER THE HARVEST RATE OF GLUCOSE, SINCE THE HARVEST 
%               RATE IS A QUANTITY GENERALLY RELATED TO THE PRODUCT OF 
%               INTEREST YIELD.

function local_sensitivity = local_sensitivity_medium_variability(  model, ...
                                                                    met_names, ...
                                                                    nominal_medium, ...
                                                                    nominal_concentration, ...
                                                                    index_data)

    F = model.F;
    n_u = length(nominal_medium); n_c = length(nominal_concentration);
    I_u = eye(n_u);
    F_u = zeros(n_c, n_u);
    F_u(index_data.decision_metabolite_indices,:) = I_u*F;

    JcQ = model.Amac*MacroKineticsJacobian(model.theta_matrix, nominal_concentration);  % Jacobian of rates with respect to concentrations
    JcM = (model.Xv*JcQ(index_data.coupled_met_indices,:) - F*eye(n_c));                % Jacobian of Mass-Balance Equation w.r.t. concentrations
    nominal_rate = ComputeQExt(nominal_concentration,model.theta_matrix,model.Amac);    
    if rank(JcM) == size(JcM,1)
        JuC = -inv(JcM)*F_u;
    else
        fprintf('Error: singular Jacobian of MBE with respect to concentrations. \n');
        local_sensitivity = NaN;
        return
    end
    JuQ = JcQ*JuC; % Jacobian of rates with respect to medium composition
    
    n_q = length(nominal_rate);
    local_sensitivity = cell(max([n_c,n_q]),1);

    local_sensitivity = evaluate_local_sensitivities( local_sensitivity,...
                                                                    JuC,...
                                                                    JuQ,...
                                                                    nominal_rate,...
                                                                    met_names,...
                                                                    model,...
                                                                    index_data);
    
end

function local_sensitivity = evaluate_local_sensitivities( local_sensitivity,...
                                                                        JuC,...
                                                                        JuQ,...
                                                                        nominal_rate,...
                                                                        met_names,...
                                                                        model,...
                                                                        index_data)
    n_c = size(JuC,1);
    n_q = size(JuQ,1);
    F = model.F;
    nabla_u_growth = JuQ(index_data.growth_index,:)';
    for k = 1 : max([n_c,n_q])
        if k <= n_c && k <= n_q
            local_sensitivity{k}.name = met_names{k};
            local_sensitivity{k}.concentrations = JuC(k,:)';
            local_sensitivity{k}.rates = JuQ(k,:)';
            local_sensitivity{k}.harvest = -nabla_u_growth*nominal_rate(k)+(F-nominal_rate(index_data.growth_index))*JuQ(k,:)';
        elseif k > n_c && k <= n_q
            local_sensitivity{k}.name = met_names{k};
            local_sensitivity{k}.concentrations = [];
            local_sensitivity{k}.rates = JuQ(k,:)';      
            local_sensitivity{k}.harvest = -nabla_u_growth*nominal_rate(k)+(F-nominal_rate(index_data.growth_index))*JuQ(k,:)';
        elseif k > n_q && k <= n_c % this will never happen I think 
            local_sensitivity{k}.name = met_names{k};
            local_sensitivity{k}.concentrations = JuC(k,:)';
            local_sensitivity{k}.rates = [];   
            local_sensitivity{k}.harvest = [];
        else
            fprintf('Error in the evaluation of the sensitivity. Index out of range. \n');
        end
    end
end
