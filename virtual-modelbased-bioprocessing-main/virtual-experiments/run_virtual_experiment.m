%% Run a virtual experiment
%   author        : Mirko Pasquini 
%
%   The function allows the user to perform a simulated experiments, once
%   the medium composition is specified (in the vector u) for all the
%   coupled metabolites. Please refer to the "Beginner's guide" for a 
%   definition of coupled metabolite.   
%
%   [c, q, w] = RUN_VIRTUAL_EXPERIMENT(u, model, c0, coupled_met_indices, display_option, runaway_res, max_iter)
%
%   @inputs:
%       u : vector of the extracellular metabolites concentrations in the 
%           medium, corresponding to the coupled metabolites.
%       model : this is a structure containing all the useful information 
%           on the kinetic model. Please refer to the "Beginner's guide" 
%           for further info.
%       c0  : vector of initial concentration, used as initialization for 
%           the solver that solves the mass-balance equation
%       coupled_met_indices : vector of indexes for the coupled
%           metabolites. The order refers to the order of metabolites in the
%           Amac matrix. Please refer to the "Beginner's guide" for a 
%           definition of coupled metabolite.
%       display_option : see optimoptions for fmincon. The default will be 
%           'none', to hide any output from the solver. Another useful 
%           option is 'iter-detailed' showing details for the iterations of the 
%           mass-balance equation solver. 
%       runaway_res : as soon as one of the solutions (obtained by starting
%           from different initial conditions around c0) has residuals
%           lower than runaway_res, the function is terminated and that
%           solution is returned (i.e. the solution is good enough). The
%           default value if the user does not input any is 1e-4. This is
%           not implemented yet.
%       max_iter : maximum number of iteration for the solver to run (if
%           left empty the default value is 50).
%
%   @outputs:
%       c : vector of concentrations at steady-state
%       q : vector of uptake-secretion rates at steady-state
%       w : vector of macroreaction rates at steady-state

function [c, q, w] = run_virtual_experiment(u, model, c0, coupled_met_indices, display_option, runaway_res, max_iter)

    if nargin < 7
        max_iter = 50;
        if nargin < 6
            runaway_res = 1e-4;
            if nargin < 5
                display_option = 'none';
            end
        end
    end

    optionsLsqNonlin = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display',display_option,'OptimalityTolerance',1e-6,'MaxIterations',max_iter,'Algorithm','levenberg-marquardt'); %  ,'Algorithm','trust-region-reflective' maybe comment out optimality tolerance, add ,'Diagnostics','on' if needed
    [c_best,resnorm_best] = lsqnonlin(@(x)MBEObj(x,u,model.F,model.Xv,model.Amac,model.theta_matrix,coupled_met_indices),c0,zeros(length(c0),1),[],optionsLsqNonlin);

    c = c_best;
    [q,w] = ComputeQExt(c,model.theta_matrix,model.Amac);
    if strcmp(display_option,'none') == false
        disp(['Mass-Balance Equation solved with 2-norm error =',' ',   num2str(norm(model.Xv*q(coupled_met_indices) - model.F*c + model.F*u),2)])
    end
end

function [obj, gradobj] = MBEObj(x,u,F,Xv,Amac,theta_matrix,coupled_met_indices)
    q = ComputeQExt(x,theta_matrix,Amac);
    q = q(coupled_met_indices);
    obj = (Xv*q - F*x + F*u);
    Jw = MacroKineticsJacobian(theta_matrix,x);
    gradobj = (Xv*Amac(coupled_met_indices,:)*Jw-F*eye(length(x)));
end
