%% Evaluate the Jacobian of the concentrations with respect to the kinetic 
%% parameters
%   author        : Mirko Pasquini
%
%  ===== NOTE : THIS SEEMS TO BE SLIGHTLY LESS EFFICIENT THAN
%  Jacobian_parameters_of_concentrations.m HOWEVER THIS IS
%  FOLLOWING A MORE ALIGNED STANDARD WITH THE REST. NEED TO CHECK WHICH 
%  FUNCTION IS IN USE. =====
%
%  This function returns the Jacobian of the extracellular metabolites
%  concentrations with respect to the kinetic parameters. This is done
%  through the implicit function theorem, using the mass-balance equation
%  as an implicit map.
%
%  JthC = EVALUATE_THETA_JACOBIAN_OF_CONCENTRATIONS(model,c,coupled_met_indices) 
%
%  @inputs: 
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide"
%           for further info.
%       c : vector of extracellular metabolite concentrations (> 0). These 
%           corresponds to coupled metabolites (Please refer to the 
%           "Beginner's guide" for a definition of coupled metabolite")
%       coupled_met_indices : vector of indexes for the coupled metabolites.
%           The order refers to the order of metabolites in the Amac matrix. 
%           Please refer to the "Beginner's guide" for a definition of 
%           coupled metabolite.
%
%   @outputs:
%       JthC : Jacobian of the extracellular metabolites concentrations, 
%           with respect to the kinetic parameters. Considering the kinetic
%           parameters vector theta (refer to the repository Wiki Page 
%           "Macroreaction Kinetics Parameters : Matrix and Vector" for 
%           a description of the structure of such vector), and the vector
%           c of extracellular (coupled) metabolite concentrations, then
%           JthC(i,j) = d c(i) / d theta(j).


function JthC = evaluate_theta_Jacobian_of_concentrations(model,c,coupled_met_indices)
    Amac = model.Amac;
    Xv = model.Xv;
    F = model.F;
    theta = model.theta_matrix;
    JthC = -inv(Xv*Amac(coupled_met_indices,:)*MacroKineticsJacobian(theta,c)-F*eye(length(c)))*...
            (Xv*Jacobian_parameters_of_rates(theta,c,Amac(coupled_met_indices,:))); % see Implicit Function Theorem
end