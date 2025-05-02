%   [EFMs,Macro_rates,q_CG_train,q_CG_predict] = Step_2_EFM_computation(data_qext,data_stoich,options_EFM,options_data)
%
%   This code performs several functons :
%   - checks if the files (.mat, .xlsx or .xls) containing the cell specific rate
%     data qext, the concentration data cext and the metabolic network exist. The directories of these files must be defined
%     in MAIN_MODELING with the following corresponding variables: directory_file_cell_specific_rate_data_qext,
%     directory_file_concentration_data_cext and directory_file_metabolic_network
%   - loads/extracts necessary informations for the modeling (such as stroichiometric matrix, the number of data, reversibility of the reactions ,etc)
%   - smoothes the data cext and qext (or not) according to the choice of the user made in MAIN_MODELING.m
%   - normalizes the qext data (or not) according to the choice of the user made in MAIN_MODELING.m
%   - strores all the necessary information for modeling in the structures data_qext, data_cext and data_stoich
%
%   @inputs:
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
%       - EFMs : matrix gathering all the computed EFMs in each of its column
%       - Macro_rates : matrix containing the computed macroscopic rates. Each row corresponds to a macroscopic reaction and the oclumns to the different experimental conditions
%       - q_CG_train : matrix contaning the output of the column generation model (without kinetics) evaluated with the training data
%       - q_CG_predict : matrix contaning the output of the column generation model (without kinetics) evaluated with the prediction data

function [EFMs,Macro_rates,q_CG_train,q_CG_predict] = Step_2_EFM_computation(data_qext,data_stoich,options_EFM,options_data)
 
clc
fprintf("Step 2: Initial test of the metabolic network...\n")
pause(0.5)

options_qp = optimoptions('quadprog','Display','none'); % 
average_data_per_condition = options_data.average_data_per_condition;

%% Load the necessary data
index_training_media = data_qext.index_training_media;
index_prediction_media = data_qext.index_prediction_media;
use_mosek = options_EFM.use_mosek;
ic = data_qext.ic;
qext_train = data_qext.qext_train';
qext_predict = data_qext.qext_predict';
Ameas = data_stoich.Ameas; 
Aint = data_stoich.Aint;
Dirrev = data_stoich.Dirrev; 
n_v = data_stoich.n_v;
n_v_irrev = data_stoich.n_v_irrev;
n_mets_int = data_stoich.n_mets_int;
[n_qext,N_train] = size(qext_train);
[~,N_predict] = size(qext_predict);

tolerance_EFM = options_EFM.tolerance_EFM;
computation_method_EFM = options_EFM.computation_method_EFM;
reduction_of_EFM = options_EFM.reduction_of_EFM;
factor_error_reduction = options_EFM.factor_error_reduction;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial data fit %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each condition/data kk = 1,...,N we we solve the followig problem
%       min ||qext(kk) - Ameas*v(kk)||^2 
%       such that Aint*v(kk) = 0            %  consequence of the pseudo-steady state assumption
%                 Dirrev*v(kk) >= 0         %  irreversibility constraint of some fluxes
% This allows to give the best data fitting peformances allowed by the
% metabolic network without any kinetic estimation.

% This step is used in order to detect a problem with the matbolic network.
% If a metabolite has a fit less than 50%, then the code stops to warn the user. The user can choose to continue by ignoring the warning.

normalization_matrix = data_qext.normalization_matrix;
Ameas = normalization_matrix*Ameas;

Aint_sp = sparse(Aint); % Convert to sparse matrix in order to increase accuracy for optimization
Ameas_sp = sparse(Ameas); % Convert to sparse matrix in order to increase accuracy for optimization
if(use_mosek == 0)
  H = Ameas'*Ameas + 1e-10*eye(n_v); % Quadratic matrix for the Matlab function quadprog
  H = (H+H')/2; % This line guarantees that the quadratic matrix is symmetric
end
Aeq_qp = Aint_sp; % Matrix used to implement the equality constraints Aint*v = 0
beq_qp = sparse(n_mets_int,1); % Vector used to implement the equality constraints Aint*v = 0
Aineq_qp = -sparse(Dirrev);  % Matrix used to implement the inequality constraints Dirrev*v >= 0 
bineq_qp = sparse(n_v_irrev,1); % Vector used to implement the equality constraints Dirrev*v >= 0 
qext_proj = zeros(n_qext,N_train); % Matrix containing the projected cell specific rate data onto the realisable set
parfor kk = 1:N_train  
  Yk = qext_train(:,kk); fk = sparse(-Ameas_sp'*Yk);
  if(use_mosek) % solve the quadratic problem with MOSEK
    v_opt_kk = quadprog_mosek_as_second_order_programming(Ameas_sp,Yk,Aineq_qp,bineq_qp,Aeq_qp,beq_qp,-Inf*ones(n_v,1),Inf*ones(n_v,1));  
  else  % solve the quadratic problem with quadprog
    v_opt_kk = quadprog(H,fk,Aineq_qp,bineq_qp,Aeq_qp,beq_qp,[],[],[],options_qp);
  end
  qext_proj(:,kk) = Ameas_sp*v_opt_kk; 
end

%% We compute the fitsand ask the user if the modeling should pursue in case of fits smaller than 50%
fit_for_each_qext = zeros(n_qext,1);
for ii = 1:n_qext
  fit_for_each_qext(ii,1) = 100*(1-norm(qext_train(ii,:) - qext_proj(ii,:),2)/norm(qext_train(ii,:) - mean(qext_train(ii,:)),2));
end
mets_meas_in_qext = data_qext.mets_meas_in_qext;
bad_fits = (fit_for_each_qext < 50);
bad_mets = mets_meas_in_qext(bad_fits);
if(isempty(bad_mets))
  fprintf("All the metabolites have a fit above 50%%.\n")
  pause(2)
  fprintf("The metabolic network seems suited for kinetics identification.\n")
  pause(2)
  fprintf("The identification will continue.\n")
  close all
else
  fprintf("These metabolites listed below have a fit less than 50%% : ")
  for ii = 1:length(bad_mets)
    if(ii == length(bad_mets))  
      fprintf("%c",bad_mets{ii});
      fprintf("\n");
    else
      fprintf("%c",bad_mets{ii});
      fprintf(", ");
    end
  end
  pause(2)
  fprintf("The metabolic network seems not suited for kinetics identification or the data are too noisy.\n")
  pause(2)
  fprintf("Press any key to continue the identification. Othwerise, stop the code (Ctrl + C).\n")
  pause
  close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Classical Column generation %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A first EFM is computed using the data
E1 = first_EFM_computation(qext_proj,Ameas,Aint,Dirrev);

switch computation_method_EFM
  % Global EFM generation  
  % All the data are considered together for the computation of all the EFMs
  case 'global' 
    clc
    text_to_print = "Step 2: EFMs computation - Global approach.\n";
    Y = qext_proj; K = Ameas; 
    [EFMs,Macro_rates,~] = EFM_column_generation(Y,K,E1,Aint,Dirrev,tolerance_EFM,use_mosek);

  % Sequential EFM generation
  % We compute an optimal subset of EFM for each data in a sequential manner
  % We compute the subset of EFM fitting the first measurement, then we try to fit the second one with the previously computed EFM and complement it if needed, etc 
  case 'sequential' 
    K = Ameas; EFMs = E1;  
    for kk = 1:N_train
      clc
      fprintf("Step 2: EFMs computation - Sequential approach.\n")    
      fprintf("Data %d/%d.\n",kk,N_train)    
      Y = qext_proj(:,kk);
      [EFMs,Macro_rates,~] = EFM_column_generation(Y,K,EFMs,Aint,Dirrev,tolerance_EFM,use_mosek);
    end

  % Unionized EFM generation
  % We compute an optimal subset of EFM for each data independently
  % All these EFMs are combined together by computing their union  
  case 'union' 
    K = Ameas; E_all = E1; 
    for kk = 1:N_train
      clc
      fprintf("Step 2: EFMs computation - Union approach.\n")    
      fprintf("Data %d/%d.\n",kk,N_train)    
      Y = qext_proj(:,kk);
      [EFMs,~,~] = EFM_column_generation(Y,K,E1,Aint,Dirrev,tolerance_EFM,use_mosek);
      E_all = combine_two_sets_of_EFMs(E_all,EFMs);    
    end
    [~,~,~,Macro_rates] = master_and_subproblem(qext_proj,K,E_all,Aint,Dirrev,use_mosek);
    EFMs = E_all;
  otherwise
  error("The variable computation_method_EFM must be a char equal to either global, sequential or union.\n")
end

% Reduction of the EFMs (if chosen by the user)
if(reduction_of_EFM) 
  Y = qext_proj;  
  [EFMs,Macro_rates,~] = reduction_EFM(Y,K,EFMs,factor_error_reduction,use_mosek);
end

negligible_wb = (sum(Macro_rates < 1e-6,2) == size(Macro_rates,2)); % Detect the possible negligible rates 
EFMs(:,negligible_wb) = []; % Remove the EFMs with negligible rates 
Macro_rates(negligible_wb,:) = []; % Remove the negligible rates
Amac = Ameas*EFMs; % Macroscopic stoichiometric matrix
Amac(abs(Amac) < 1e-8) = 0; % Macroscopic the negligible stoichiometric coefficients

if(N_predict > 0) % Compute the rates for the prediction media
  [~,~,~,Macro_rates_predict] = master_and_subproblem(qext_predict,K,EFMs,Aint,Dirrev,use_mosek);
else
  Macro_rates_predict = [];  
end

% If averaging has been selected, this part of the code just simply
% duplicates the value obtained for each condition for the corresponding
% measurement time instants
if(average_data_per_condition == 1) 
  W_CG_train = []; W_CG_predict = [];
  kl_train = 1;
  for k = index_training_media
    ind_k = (ic == k); nb_data_in_condition_k = sum(ind_k);
    for l = 1:nb_data_in_condition_k
      W_CG_train = [W_CG_train,Macro_rates(:,kl_train)];
    end
    kl_train = kl_train + 1;
  end
  if(~isempty(index_prediction_media))
    kl_predict = 1;
    for k = index_prediction_media
      ind_k = (ic == k); nb_data_in_condition_k = sum(ind_k);
      for l = 1:nb_data_in_condition_k
        W_CG_predict = [W_CG_predict,Macro_rates_predict(:,kl_predict)];
      end
      kl_predict = kl_predict + 1;
    end
  end
  w_CG_predict =  W_CG_predict; w_CG_train =  W_CG_train;
else
  w_CG_predict = Macro_rates_predict;  
  w_CG_train = Macro_rates;  
end

if(N_predict > 0)
  q_CG_predict = Amac*w_CG_predict; % Cell specific rate model without kinetics evaluated with the prediction data
else
  q_CG_predict = [];
end

q_CG_train = Amac*w_CG_train; % Cell specific rate model without kinetics evaluated with the training data

end

function [E1] = first_EFM_computation(Y,K,Aint,Dirrev)
   % This code computes a first EFM from the matrix data in Y such that the least squares error Y-K*E*w is reduced  
   % It solves the following linear programming where e is the EFM vector (variable)
   %             min -Y'*K*e
   %   such that Aint*e = 0
   %             Dirrev*e >= 0
   %             ||e||_1 = 1

   options_lin = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
   options_lin.Preprocess = 'none';
   f0_all = -K'*Y; % vector of the linear objective function
   fo = sparse(sum(f0_all,2)); % conversion of the objective vector into a sparse vector for accuracy improvement 
   [n_v_irrev,n_v] = size(Dirrev); 
   n_mets_int = size(Aint,1);
   
   A_lin = [-sparse(eye(n_v)),sparse(eye(n_v));
             sparse(eye(n_v)),sparse(eye(n_v));
             Dirrev,sparse(n_v_irrev,n_v)]; % matrix which implements all the linear inequalities in A_lin*x >= b_lin
   b_lin = [sparse(2*n_v+n_v_irrev,1)];
   A_eq = [Aint,sparse(n_mets_int,n_v);sparse(1,n_v),ones(1,n_v)];
   b_eq = [sparse(n_mets_int,1);1];
   f = [fo;zeros(n_v,1)];
   e_opt = linprog(f,-A_lin,-b_lin,A_eq,b_eq,[],[],options_lin); % computation of the EFM
   [e_opt2] = improve_EFM_computation(e_opt,A_lin,b_lin,A_eq,b_eq,f); % improvement of the accuracy of the computed EFM
   while(sum((e_opt2) > 0) < sum((e_opt) > 0)) % while there are always zero entries in the EFM, continue improving the accuracy by removines the corresponding variables
      e_opt = e_opt2; 
      [e_opt2] = improve_EFM_computation(e_opt,A_lin,b_lin,A_eq,b_eq,f);
   end
   E1 = e_opt2(1:n_v);
   
   if(abs(norm(E1,1) - 1) > 1e-3)
     error("No EFMs could befound for initialization.\n")
   end

end

function [e_opt_new] = improve_EFM_computation(e_opt,A_lin,b_lin,A_eq,b_eq,f)

% The purpose of this code is to compute again an EFM where all the entries which are equal to 0 are removed
% The idea is to iteratively reduce the EFM computation complexity (reduce the number of variables)
% It will increase the accuracy of the new optimum
%
%   @inputs:  
%       - e_opt : previously computed EFM for which we want to improve the accuracy (before removing the zero entries)
%       - A_lin : matrix which implements the linear constraints A_lin*x >= b_lin where x is the variable vector
%       - b_lin : vector which implements the linear constraints A_lin*x >= b_lin where x is the variable vector
%       - A_eq : matrix which implements the linear constraints A_eq*x = b_eq where x is the variable vector
%       - f : vector which implements the linear objective function f'*x where x is the variable vector  
% 
%   @outputs:  
%       - e_opt_new : new EFM matrix for which the accuracy has possibly been improved

  options_lin = optimoptions('linprog','Display','none');
  options_lin.Preprocess = 'none';
  ne = length(e_opt); e_opt_new = zeros(ne,1);
  ind_e_nonzeros = nonzeros(diag(abs(e_opt) > 1e-9)*(1:ne)'); % compute the indices of the entries of e_opt which are non-zero
  A_lin_red = A_lin(:,ind_e_nonzeros); % keep the columns of A_lin corresponding to the variables which are non-zeros
  A_eq_red = A_eq(:,ind_e_nonzeros); % keep the columns of A_eq corresponding to the variables which are non-zeros
  f_red = f(ind_e_nonzeros); % keep the entries of f corresponding to the variables which are non-zeros
  e_opt_red = linprog(f_red,-A_lin_red,-b_lin,A_eq_red,b_eq,[],[],options_lin); % compute the non-zeros entries of the EFM e_opt
  e_opt_new(ind_e_nonzeros) = e_opt_red; % form the final EFM with the zero entries included

end

function [Eb,Wb,cost_MP] = reduction_EFM(Y,K,Eb,factor_error_reduction,use_mosek)

  % This code proposes to reduce the number of EFMs computed after the
  % column generation algorithm.
  % The idea of the algorithm is as follows
  % - Step 0: we compute the minimal cost with the whole subset of EFMs
  % - Step 1: for each EFM in Eb, we remove it from the subset and we
  %           compute the least-squares error fit by solving the
  %           masterproblem. 
  % - Step 2: among all these errors, we compute the minimal one.
  % - Step 3: if this minimal cost is smaller than the threshold cost, we remove
  %           the corresponding EFM and we re-iteratre Steps 1, 2 and 3. If not, 
  %           the reduction is over. 

  % The threshold cost is computed as follows
  %     cost_MP_threshold = (100 + factor_error_reduction)/100*cost_MP
  % where cost_MP is the error cost obtained in Step 0 with the whole
  % subset of EFMs and factor_error_reduction is the pourcentage of cost
  % degradation allowed by the user for EFM reduction. For instance, if the
  % user sets factor_error_reduction = 45, then this code will remove as
  % many EFMs as possible as long as the initial error cost has not
  % increased more than 45% from its original value (with the whole subset of
  % EFMs)

  options_qp = optimoptions('quadprog','Display','none');   

  fprintf("Reduction of the EFMs.\n")   
  
  % Step 0: We first compute the least-squares error fit of the whole subset of EFM
  % with respect to the data
  % For this purpose, we compute the Macro-rates by solving a quadratic
  % cost (masterproblem cost)
  N = size(Y,2);
  n_Eb = size(Eb,2);
  L = K*Eb; 
  H = L'*L + 1e-15*eye(n_Eb); 
  H = (H+H')/2;
  Wb = zeros(n_Eb,N); 
  cost_MP = 0; 
  parfor kk = 1:N
    Yk = Y(:,kk);  
    fk = -L'*Yk;  
    if(use_mosek)
      x = quadprog_mosek_as_second_order_programming(L,Yk,[],[],[],[],zeros(n_Eb,1),Inf*ones(n_Eb,1));
    else
      x = quadprog(H,fk,[],[],[],[],zeros(n_Eb,1),[],[],options_qp);
    end
    x(x < 0) = 0;
    Wb(:,kk) = x;
    cost_MP = cost_MP + norm(Y(:,kk)-L*Wb(:,kk),2)^2;
  end

  cost_MP_threshold = (100 + factor_error_reduction)/100*cost_MP;  % cost_MP_threshold is the value of the least-squares error fit
  reduction_completed = 0;
  index_EFM_removed = [];
  while(reduction_completed == 0)    
    % We remove the EFM with index index_EFM_removed which leads to a minimal cost lower than
    % the threshold cost cost_MP_threshold. At the first iteration, we do not remove any EFM.

    Eb(:,index_EFM_removed) = [];
    Wb(index_EFM_removed,:) = [];
    % W look for negligible Wb and remove the corresponding EFMs  
    % negligible_wb = (sum(Wb < 1e-6,2) == N);
    % Eb(:,negligible_wb) = []; 

    % We start the reduction process here
    n_Eb = size(Eb,2);
    all_cost_MP_red = zeros(n_Eb,1);
    clc
    fprintf("Step 2: EFMs computation - Reduction of the EFMs.\n")   
    fprintf("Number of remaining EFMs: %d.\n",n_Eb)
    for j = 1:n_Eb % Step 1: for each EFM in Eb, we compute the error obtained by removing one EFM at a time. 
      Eb_red = Eb; Eb_red(:,j) = []; % New subset of EFM with the j-th EFM removed from Eb
      n_Eb_red = size(Eb_red,2); 
      L = K*Eb_red; 
      if(~use_mosek)
        H = L'*L + 1e-15*eye(n_Eb_red); 
        H = (H+H')/2;
      end
      Wb_red = zeros(n_Eb_red,N); 
      cost_MP_red = 0;
      parfor kk = 1:N  % Masterproblem with one EFM removed
        Yk = Y(:,kk);  
        fk = -L'*Yk;  
        if(use_mosek)
          x = quadprog_mosek_as_second_order_programming(L,Yk,[],[],[],[],zeros(n_Eb_red,1),Inf*ones(n_Eb_red,1));  
        else
          x = quadprog(H,fk,[],[],[],[],zeros(n_Eb_red,1),[],[],options_qp);
        end
        x(x < 0) = 0;
        Wb_red(:,kk) = x;
        cost_MP_red = cost_MP_red + norm(Y(:,kk)-L*Wb_red(:,kk),2)^2;
      end
      all_cost_MP_red(j,1) = cost_MP_red; % Error obtained with the j-th EFM removed
    end
    [cost_MP_red_min,index_EFM_removed] = min(all_cost_MP_red); %  Step 2: among all these errors, we compute the minimal one. We also compute the corresponding index of the removed EFM whi
    if(cost_MP_red_min > cost_MP_threshold) %  Step 3: we check if the minimal cost cost_MP_red_min is above the threshold cost
      reduction_completed = 1; % If yes, we stop the reduction, If no: we continue
    end
  end

  

  % We compute the macroscopic rates and error obtained with the reduced
  % subset of EFM by solving once again the masterproblem
  n_Eb = size(Eb,2);
  L = K*Eb; 
  if(~use_mosek)
    H = L'*L + 1e-15*eye(n_Eb); 
    H = (H+H')/2;
  end
  Wb = zeros(n_Eb,N); 
  cost_MP = 0; 
  parfor kk = 1:N
    Yk = Y(:,kk);  
    fk = -L'*Yk;  
    if(use_mosek)
      x = quadprog_mosek_as_second_order_programming(L,Yk,[],[],[],[],zeros(n_Eb,1),Inf*ones(n_Eb,1));  
    else
      x = quadprog(H,fk,[],[],[],[],zeros(n_Eb,1),[],[],options_qp);
    end
    x(x < 0) = 0;
    Wb(:,kk) = x;
    cost_MP = cost_MP + norm(Y(:,kk)-L*Wb(:,kk),2)^2;
  end

  negligible_wb = (sum(Wb < 1e-6,2) == N);
  Eb(:,negligible_wb) = [];
  Wb(negligible_wb,:) = [];
  n_Eb = size(Eb,2);
  clc
  fprintf("Step 2: EFMs computation - Reduction of the EFMs.\n")   
  fprintf("Number of remaining EFMs: %d.\n",n_Eb)
  pause(2)
end

function [E_union] = combine_two_sets_of_EFMs(E1,E2)

%% This function computes the union between two sets of EFMs.
%
%   @inputs:  
%       - E1 : first EFM matrix for which the union will be computed from with E2
%       - E2 : second EFM matrix for which the union will be computed from with E1
%
%   @outputs:  
%       - E_union : EFM matrix which is the union of E1 and E2

E_union = E1;
n1 = size(E1,2); % number of EFMs in E1
n2 = size(E2,2); % number of EFMs in E2

for i = 1:n2 % For each EFM in E2, we check if it exists in E1
  already_in_set = 0;
  for j = 1:n1 % We try each EFM in E1 and check if it is teh same as the i-th EFM in E2
    if(abs(E1(:,j)-E2(:,i)) < 1e-6)
      already_in_set = 1; % If we find the same EFM, we stop the search among all the EFMs in E2
      break
    end
  end
  if(already_in_set == 1) % If the same EFM has been found, we do not include it in E_union
  else
    E_union = [E_union,E2(:,i)];
  end
end

end
